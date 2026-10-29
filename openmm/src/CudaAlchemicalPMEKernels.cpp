#include "AlchemicalPMEKernels.h"

#include "openmm/KernelFactory.h"
#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "AlchemicalPMEUtilities.h"
#include "openmm/internal/ContextImpl.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#ifdef OPENMM_COLVARS_ENABLE_CUDA_KERNELS

#if __has_include("CudaArray.h")
#include "CudaArray.h"
#include "CudaContext.h"
#include "CudaPlatform.h"
#else
#include "openmm/cuda/CudaArray.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/cuda/CudaPlatform.h"
#endif
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeForceInfo.h"
#include "openmm/common/ComputeParameterInfo.h"
#include "openmm/common/FFT3D.h"
#include "openmm/common/NonbondedUtilities.h"
#endif

using namespace OpenMM;

namespace {

#ifdef OPENMM_COLVARS_ENABLE_CUDA_KERNELS

const char* AlchemicalPMECudaSourcePart1 = R"(
struct ParticleParam {
    real charge0;
    real chargeDelta;
    real sigma0;
    real sigmaDelta;
    real epsilon0;
    real epsilonDelta;
};

DEVICE real3 make_real3(real v) {
    return make_real3(v, v, v);
}

DEVICE real2 make_real2(real v) {
    return make_real2(v, v);
}

DEVICE real3 operator+(real3 a, real3 b) {
    return make_real3(a.x+b.x, a.y+b.y, a.z+b.z);
}

DEVICE real3 operator-(real3 a, real3 b) {
    return make_real3(a.x-b.x, a.y-b.y, a.z-b.z);
}

DEVICE real3 operator-(real3 a) {
    return make_real3(-a.x, -a.y, -a.z);
}

DEVICE real3 operator*(real3 a, real3 b) {
    return make_real3(a.x*b.x, a.y*b.y, a.z*b.z);
}

DEVICE real3 operator*(real a, real3 b) {
    return make_real3(a*b.x, a*b.y, a*b.z);
}

DEVICE real3 operator*(real3 a, real b) {
    return make_real3(a.x*b, a.y*b, a.z*b);
}

DEVICE double3 makeDouble3(double x, double y, double z) {
    double3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

DEVICE double dot3(double3 a, double3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

DEVICE double3 cross3(double3 a, double3 b) {
    return makeDouble3(a.y*b.z - a.z*b.y,
                       a.z*b.x - a.x*b.z,
                       a.x*b.y - a.y*b.x);
}

DEVICE double3 scale3(double3 a, double s) {
    return makeDouble3(a.x*s, a.y*s, a.z*s);
}

DEVICE double3 add3(double3 a, double3 b) {
    return makeDouble3(a.x+b.x, a.y+b.y, a.z+b.z);
}

DEVICE double3 sub3(double3 a, double3 b) {
    return makeDouble3(a.x-b.x, a.y-b.y, a.z-b.z);
}

DEVICE real clamp01(real value) {
    return fmin((real) 1, fmax((real) 0, value));
}

DEVICE void interpolateElectrostaticEndpoint(real value0, real valueDelta, real lambda,
                                             real electrostaticLambdaStart,
                                             real* value, real* dValue) {
    if (electrostaticLambdaStart <= (real) 0) {
        *value = value0 + lambda*valueDelta;
        *dValue = valueDelta;
        return;
    }
    const real value1 = value0 + valueDelta;
    const real valueTol = (real) 1.0e-12;
    const bool grow = fabs(value0) <= valueTol && fabs(value1) > valueTol;
    const bool vanish = fabs(value1) <= valueTol && fabs(value0) > valueTol;
    if (grow) {
        const real denom = (real) 1 - electrostaticLambdaStart;
        if (lambda <= electrostaticLambdaStart) {
            *value = (real) 0;
            *dValue = (real) 0;
        }
        else {
            const real schedule = clamp01((lambda - electrostaticLambdaStart)/denom);
            *value = value1*schedule;
            *dValue = schedule < (real) 1 ? value1/denom : (real) 0;
        }
        return;
    }
    if (vanish) {
        const real end = (real) 1 - electrostaticLambdaStart;
        if (lambda >= end) {
            *value = (real) 0;
            *dValue = (real) 0;
        }
        else {
            const real schedule = clamp01((real) 1 - lambda/end);
            *value = value0*schedule;
            *dValue = schedule > (real) 0 ? -value0/end : (real) 0;
        }
        return;
    }
    *value = value0 + lambda*valueDelta;
    *dValue = valueDelta;
}

DEVICE void saveForce(int atom, double3 force, GLOBAL mm_ulong* RESTRICT forceBuffers) {
    ATOMIC_ADD(&forceBuffers[atom], (mm_ulong) realToFixedPoint((real) force.x));
    ATOMIC_ADD(&forceBuffers[atom+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint((real) force.y));
    ATOMIC_ADD(&forceBuffers[atom+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint((real) force.z));
}

)"; 

const char* AlchemicalPMECudaSourcePart2 = R"(
KERNEL void findAlchemicalPMEAtomGridIndex(
        GLOBAL const real4* RESTRICT posq,
        GLOBAL int2* RESTRICT pmeAtomGridIndex,
        real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ) {
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % GRID_SIZE_X,
                                   ((int) t.y) % GRID_SIZE_Y,
                                   ((int) t.z) % GRID_SIZE_Z);
        pmeAtomGridIndex[atom] = make_int2(atom, gridIndex.x*GRID_SIZE_Y*GRID_SIZE_Z+gridIndex.y*GRID_SIZE_Z+gridIndex.z);
    }
}

KERNEL void spreadAlchemicalPMECharge(
        GLOBAL const real4* RESTRICT posq,
        GLOBAL real* RESTRICT pmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ,
        GLOBAL const int2* RESTRICT pmeAtomGridIndex,
        GLOBAL const int* RESTRICT atomIndex,
        GLOBAL const real* RESTRICT charges) {
    real3 data[PME_ORDER];
    const real scale = RECIP((real) (PME_ORDER-1));
    for (int i = GLOBAL_ID; i < NUM_ATOMS*PME_ORDER; i += GLOBAL_SIZE) {
        int atom = pmeAtomGridIndex[i/PME_ORDER].x;
        int originalAtom = atomIndex[atom];
        real4 pos = posq[atom];
        const real charge = (real) charges[originalAtom]*EPSILON_FACTOR;
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % GRID_SIZE_X,
                                   ((int) t.y) % GRID_SIZE_Y,
                                   ((int) t.z) % GRID_SIZE_Z);
        if (charge == 0)
            continue;
        real3 dr = make_real3(t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z);
        data[PME_ORDER-1] = make_real3(0);
        data[1] = dr;
        data[0] = make_real3(1)-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP((real) (j-1));
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+make_real3(k))*data[j-k-2] + (make_real3(j-k)-dr)*data[j-k-1]);
            data[0] = div*(make_real3(1)-dr)*data[0];
        }
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+make_real3(j))*data[PME_ORDER-j-2] + (make_real3(PME_ORDER-j)-dr)*data[PME_ORDER-j-1]);
        data[0] = scale*(make_real3(1)-dr)*data[0];

        int iz = i%PME_ORDER;
        int zindex = gridIndex.z+iz;
        zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
        real dz = 0;
        for (int j = 0; j < PME_ORDER; j++)
            dz = j == iz ? data[j].z : dz;
        dz *= charge;
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = gridIndex.x+ix;
            xbase -= (xbase >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            xbase = xbase*GRID_SIZE_Y*GRID_SIZE_Z;
            real dzdx = dz*data[ix].x;
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndex.y+iy;
                ybase -= (ybase >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                ybase = ybase*GRID_SIZE_Z;
                int index = xbase + ybase + zindex;
                real add = dzdx*data[iy].y;
                if (fabs(add) > 2.3e-10f)
                    ATOMIC_ADD(&pmeGrid[index], add);
            }
        }
    }
}

)";

const char* AlchemicalPMECudaSourcePart3 = R"(
KERNEL void convolveAlchemicalPMEGrid(
        GLOBAL real2* RESTRICT pmeGrid,
        GLOBAL const real* RESTRICT pmeBsplineModuliX,
        GLOBAL const real* RESTRICT pmeBsplineModuliY,
        GLOBAL const real* RESTRICT pmeBsplineModuliZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ) {
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*(GRID_SIZE_Z/2+1);
    const real recipScaleFactor = RECIP(M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
    for (int index = GLOBAL_ID; index < gridSize; index += GLOBAL_SIZE) {
        int kx = index/(GRID_SIZE_Y*(GRID_SIZE_Z/2+1));
        int remainder = index-kx*GRID_SIZE_Y*(GRID_SIZE_Z/2+1);
        int ky = remainder/(GRID_SIZE_Z/2+1);
        int kz = remainder-ky*(GRID_SIZE_Z/2+1);
        int mx = (kx < (GRID_SIZE_X+1)/2) ? kx : (kx-GRID_SIZE_X);
        int my = (ky < (GRID_SIZE_Y+1)/2) ? ky : (ky-GRID_SIZE_Y);
        int mz = (kz < (GRID_SIZE_Z+1)/2) ? kz : (kz-GRID_SIZE_Z);
        real mhx = mx*recipBoxVecX.x;
        real mhy = mx*recipBoxVecY.x+my*recipBoxVecY.y;
        real mhz = mx*recipBoxVecZ.x+my*recipBoxVecZ.y+mz*recipBoxVecZ.z;
        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        real bx = pmeBsplineModuliX[kx];
        real by = pmeBsplineModuliY[ky];
        real bz = pmeBsplineModuliZ[kz];
        real2 grid = pmeGrid[index];
        if (kx != 0 || ky != 0 || kz != 0) {
            real denom = m2*bx*by*bz;
            real eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;
            pmeGrid[index] = make_real2(grid.x*eterm, grid.y*eterm);
        }
        else {
            pmeGrid[index] = make_real2(0);
        }
    }
}

KERNEL void evaluateAlchemicalPMEReciprocalEnergy(
        GLOBAL const real2* RESTRICT pmeGrid,
        GLOBAL const real* RESTRICT pmeBsplineModuliX,
        GLOBAL const real* RESTRICT pmeBsplineModuliY,
        GLOBAL const real* RESTRICT pmeBsplineModuliZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ,
        GLOBAL mixed* RESTRICT energyBuffer) {
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    const real recipScaleFactor = RECIP(M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
    mixed energy = 0;
    for (int index = GLOBAL_ID; index < gridSize; index += GLOBAL_SIZE) {
        int kx = index/(GRID_SIZE_Y*GRID_SIZE_Z);
        int remainder = index-kx*GRID_SIZE_Y*GRID_SIZE_Z;
        int ky = remainder/GRID_SIZE_Z;
        int kz = remainder-ky*GRID_SIZE_Z;
        int mx = (kx < (GRID_SIZE_X+1)/2) ? kx : (kx-GRID_SIZE_X);
        int my = (ky < (GRID_SIZE_Y+1)/2) ? ky : (ky-GRID_SIZE_Y);
        int mz = (kz < (GRID_SIZE_Z+1)/2) ? kz : (kz-GRID_SIZE_Z);
        real mhx = mx*recipBoxVecX.x;
        real mhy = mx*recipBoxVecY.x+my*recipBoxVecY.y;
        real mhz = mx*recipBoxVecZ.x+my*recipBoxVecZ.y+mz*recipBoxVecZ.z;
        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        real bx = pmeBsplineModuliX[kx];
        real by = pmeBsplineModuliY[ky];
        real bz = pmeBsplineModuliZ[kz];
        real denom = m2*bx*by*bz;
        real eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;
        if (kz >= (GRID_SIZE_Z/2+1)) {
            kx = ((kx == 0) ? kx : GRID_SIZE_X-kx);
            ky = ((ky == 0) ? ky : GRID_SIZE_Y-ky);
            kz = GRID_SIZE_Z-kz;
        }
        int halfIndex = kz + ky*(GRID_SIZE_Z/2+1)+kx*(GRID_SIZE_Y*(GRID_SIZE_Z/2+1));
        real2 grid = pmeGrid[halfIndex];
        if (kx != 0 || ky != 0 || kz != 0)
            energy += (mixed) eterm*((mixed) grid.x*grid.x + (mixed) grid.y*grid.y);
    }
    energyBuffer[GLOBAL_ID] = (mixed) 0.5*energy;
}

KERNEL void computeAlchemicalPMEChargesAndScalarTerms(
        GLOBAL const ParticleParam* RESTRICT particleParams,
        real lambda,
        real electrostaticLambdaStart, real alpha,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        GLOBAL real* RESTRICT charges,
        GLOBAL real* RESTRICT dCharges,
        GLOBAL mixed* RESTRICT scalarBuffer) {
    extern __shared__ real tempBuffer[];
    real* selfEnergyBuffer = tempBuffer;
    real* selfDerivativeBuffer = tempBuffer + blockDim.x;
    real* totalChargeBuffer = tempBuffer + 2*blockDim.x;
    real* totalDChargeBuffer = tempBuffer + 3*blockDim.x;
    const int thread = threadIdx.x;
    real selfEnergy = 0;
    real selfDerivative = 0;
    real totalCharge = 0;
    real totalDCharge = 0;
    for (int atom = thread; atom < NUM_ATOMS; atom += blockDim.x) {
        const ParticleParam params = particleParams[atom];
        real charge = 0;
        real dCharge = 0;
        interpolateElectrostaticEndpoint(params.charge0, params.chargeDelta, lambda, electrostaticLambdaStart,
                                         &charge, &dCharge);
        charges[atom] = charge;
        dCharges[atom] = dCharge;
        selfEnergy += charge*charge;
        selfDerivative += (real) 2*charge*dCharge;
        totalCharge += charge;
        totalDCharge += dCharge;
    }
    selfEnergyBuffer[thread] = selfEnergy;
    selfDerivativeBuffer[thread] = selfDerivative;
    totalChargeBuffer[thread] = totalCharge;
    totalDChargeBuffer[thread] = totalDCharge;
    __syncthreads();
    for (int offset = blockDim.x/2; offset > 0; offset >>= 1) {
        if (thread < offset) {
            selfEnergyBuffer[thread] += selfEnergyBuffer[thread+offset];
            selfDerivativeBuffer[thread] += selfDerivativeBuffer[thread+offset];
            totalChargeBuffer[thread] += totalChargeBuffer[thread+offset];
            totalDChargeBuffer[thread] += totalDChargeBuffer[thread+offset];
        }
        __syncthreads();
    }
    if (thread == 0) {
        double3 a = makeDouble3((double) periodicBoxVecX.x, (double) periodicBoxVecX.y, (double) periodicBoxVecX.z);
        double3 b = makeDouble3((double) periodicBoxVecY.x, (double) periodicBoxVecY.y, (double) periodicBoxVecY.z);
        double3 c = makeDouble3((double) periodicBoxVecZ.x, (double) periodicBoxVecZ.y, (double) periodicBoxVecZ.z);
        const double volume = dot3(a, cross3(b, c));
        const mixed chargeSum = (mixed) totalChargeBuffer[0];
        const mixed dChargeSum = (mixed) totalDChargeBuffer[0];
        scalarBuffer[0] = -(mixed) 138.935456*alpha*selfEnergyBuffer[0]/(mixed) 1.77245385090551602729
                          -(mixed) 138.935456*M_PI*chargeSum*chargeSum/((mixed) 2*volume*alpha*alpha);
        scalarBuffer[1] = -(mixed) 138.935456*alpha*selfDerivativeBuffer[0]/(mixed) 1.77245385090551602729
                          -(mixed) 138.935456*M_PI*chargeSum*dChargeSum/(volume*alpha*alpha);
    }
}

)";

const char* AlchemicalPMECudaSourcePart4 = R"(
KERNEL void interpolateAlchemicalPMEForceAndDerivative(
        GLOBAL const real4* RESTRICT posq,
        GLOBAL mm_ulong* RESTRICT forceBuffers,
        GLOBAL const real* RESTRICT pmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ,
        GLOBAL const int2* RESTRICT pmeAtomGridIndex,
        GLOBAL const int* RESTRICT atomIndex,
        GLOBAL const real* RESTRICT charges,
        int includeForces,
        GLOBAL real* RESTRICT chargeDerivatives) {
    real3 data[PME_ORDER];
    real3 ddata[PME_ORDER];
    const real scale = RECIP((real) (PME_ORDER-1));
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        int atom = pmeAtomGridIndex[i].x;
        int originalAtom = atomIndex[atom];
        real3 force = make_real3(0);
        real derivative = 0;
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % GRID_SIZE_X,
                                   ((int) t.y) % GRID_SIZE_Y,
                                   ((int) t.z) % GRID_SIZE_Z);
        real q = (real) charges[originalAtom]*EPSILON_FACTOR;

        real3 dr = make_real3(t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z);
        data[PME_ORDER-1] = make_real3(0);
        data[1] = dr;
        data[0] = make_real3(1)-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP((real) (j-1));
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+make_real3(k))*data[j-k-2] + (make_real3(j-k)-dr)*data[j-k-1]);
            data[0] = div*(make_real3(1)-dr)*data[0];
        }
        ddata[0] = -data[0];
        for (int j = 1; j < PME_ORDER; j++)
            ddata[j] = data[j-1]-data[j];
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+make_real3(j))*data[PME_ORDER-j-2] + (make_real3(PME_ORDER-j)-dr)*data[PME_ORDER-j-1]);
        data[0] = scale*(make_real3(1)-dr)*data[0];

        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = gridIndex.x+ix;
            xbase -= (xbase >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            xbase = xbase*GRID_SIZE_Y*GRID_SIZE_Z;
            real dx = data[ix].x;
            real ddx = ddata[ix].x;
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndex.y+iy;
                ybase -= (ybase >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                ybase = xbase + ybase*GRID_SIZE_Z;
                real dy = data[iy].y;
                real ddy = ddata[iy].y;
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    real gridvalue = pmeGrid[ybase + zindex];
                    force.x += ddx*dy*data[iz].z*gridvalue;
                    force.y += dx*ddy*data[iz].z*gridvalue;
                    force.z += dx*dy*ddata[iz].z*gridvalue;
                    derivative += dx*dy*data[iz].z*gridvalue;
                }
            }
        }
        real forceX = -q*(force.x*GRID_SIZE_X*recipBoxVecX.x);
        real forceY = -q*(force.x*GRID_SIZE_X*recipBoxVecY.x+force.y*GRID_SIZE_Y*recipBoxVecY.y);
        real forceZ = -q*(force.x*GRID_SIZE_X*recipBoxVecZ.x+force.y*GRID_SIZE_Y*recipBoxVecZ.y+force.z*GRID_SIZE_Z*recipBoxVecZ.z);
        if (includeForces) {
            ATOMIC_ADD(&forceBuffers[atom], (mm_ulong) realToFixedPoint(forceX));
            ATOMIC_ADD(&forceBuffers[atom+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(forceY));
            ATOMIC_ADD(&forceBuffers[atom+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(forceZ));
        }
        chargeDerivatives[originalAtom] = derivative*EPSILON_FACTOR;
    }
}

KERNEL void reduceAlchemicalPMEReciprocalDerivative(
        GLOBAL const real* RESTRICT chargeDerivatives,
        GLOBAL const real* RESTRICT dCharges,
        GLOBAL mixed* RESTRICT scalarBuffer) {
    extern __shared__ mixed derivTempBuffer[];
    const int thread = threadIdx.x;
    mixed derivative = 0;
    for (int atom = thread; atom < NUM_ATOMS; atom += blockDim.x)
        derivative += (mixed) chargeDerivatives[atom]*dCharges[atom];
    derivTempBuffer[thread] = derivative;
    __syncthreads();
    for (int offset = blockDim.x/2; offset > 0; offset >>= 1) {
        if (thread < offset)
            derivTempBuffer[thread] += derivTempBuffer[thread+offset];
        __syncthreads();
    }
    if (thread == 0)
        scalarBuffer[2] = derivTempBuffer[0];
}

KERNEL void finishAlchemicalPMEDeviceBuffers(
        GLOBAL const mixed* RESTRICT pmeEnergyBuffer,
        GLOBAL const mixed* RESTRICT scalarBuffer,
        GLOBAL mixed* RESTRICT energyBuffer,
        GLOBAL mixed* RESTRICT energyParamDerivs,
        int pmeWorkSize,
        int energyBufferSize,
        int numDerivatives,
        int lambdaDerivativeIndex,
        int includeEnergy) {
    const int index = GLOBAL_ID;
    if (includeEnergy) {
        mixed energy = 0;
        for (int i = index; i < pmeWorkSize; i += energyBufferSize)
            energy += pmeEnergyBuffer[i];
        energyBuffer[index] += energy;
    }
    if (index == 0) {
        if (includeEnergy)
            energyBuffer[0] += scalarBuffer[0];
        if (numDerivatives > 0 && lambdaDerivativeIndex >= 0)
            energyParamDerivs[lambdaDerivativeIndex] += scalarBuffer[1] + scalarBuffer[2];
    }
}
)";

class AlchemicalPMEForceInfo final : public ComputeForceInfo {
public:
    explicit AlchemicalPMEForceInfo(const AlchemicalPMEForce& force) : force_(force) {
    }

    bool areParticlesIdentical(int particle1, int particle2) override {
        double charge01, charge11, sigma01, sigma11, epsilon01, epsilon11;
        double charge02, charge12, sigma02, sigma12, epsilon02, epsilon12;
        force_.getParticleParameters(particle1, charge01, charge11, sigma01, sigma11, epsilon01, epsilon11);
        force_.getParticleParameters(particle2, charge02, charge12, sigma02, sigma12, epsilon02, epsilon12);
        return charge01 == charge02 && charge11 == charge12 &&
               sigma01 == sigma02 && sigma11 == sigma12 &&
               epsilon01 == epsilon02 && epsilon11 == epsilon12;
    }

private:
    const AlchemicalPMEForce& force_;
};

class AlchemicalPMEDerivativeParameterMirror final : public ComputeContext::ForcePostComputation {
public:
    AlchemicalPMEDerivativeParameterMirror(ContextImpl& context,
                                           const AlchemicalPMEForce& force)
        : context_(context),
          lambdaParameterName_(force.getLambdaParameterName()),
          derivativeParameterName_(force.getLambdaDerivativeParameterName()),
          forceGroupFlag_(std::uint32_t{1} << static_cast<unsigned>(force.getForceGroup()))
    {
    }

    double computeForceAndEnergy(bool includeForces,
                                 bool includeEnergy,
                                 int groups) override
    {
        (void) includeForces;
        (void) includeEnergy;
        if ((static_cast<std::uint32_t>(groups) & forceGroupFlag_) == 0) {
            return 0.0;
        }

        auto& parameters = const_cast<std::map<std::string, double>&>(context_.getParameters());
        auto derivativeParameter = parameters.find(derivativeParameterName_);
        if (derivativeParameter == parameters.end()) {
            return 0.0;
        }

        std::map<std::string, double> derivatives;
        context_.getEnergyParameterDerivatives(derivatives);
        const auto derivative = derivatives.find(lambdaParameterName_);
        if (derivative != derivatives.end()) {
            derivativeParameter->second = derivative->second;
        }
        return 0.0;
    }

private:
    ContextImpl& context_;
    std::string lambdaParameterName_;
    std::string derivativeParameterName_;
    std::uint32_t forceGroupFlag_;
};

class CudaCalcAlchemicalPMEForceKernel final : public CalcAlchemicalPMEForceKernel {
public:
    CudaCalcAlchemicalPMEForceKernel(std::string name, const Platform& platform)
        : CalcAlchemicalPMEForceKernel(name, platform),
          owner_(nullptr),
          cu_(nullptr),
          pmeGridIndexKernel_(nullptr),
          pmeSpreadChargeKernel_(nullptr),
          pmeConvolutionKernel_(nullptr),
          pmeEnergyKernel_(nullptr),
          pmeInterpolateKernel_(nullptr),
          pmeChargeAndScalarKernel_(nullptr),
          pmeDerivativeReductionKernel_(nullptr),
          pmeFinishKernel_(nullptr),
          numAtoms_(0),
          numExceptions_(0),
          lambdaDerivativeIndex_(-1),
          lambdaGlobalParamIndex_(-1),
          threadBlockSize_(128),
          pmeGridSizeX_(0),
          pmeGridSizeY_(0),
          pmeGridSizeZ_(0),
          pmeGridElements_(0),
          pmeWorkSize_(0) {
    }

    void initialize(ContextImpl& context, const AlchemicalPMEForce& force) override {
        owner_ = &force;
        numAtoms_ = force.getNumParticles();
        numExceptions_ = force.getNumExceptions();
        if (numAtoms_ != context.getSystem().getNumParticles()) {
            throw OpenMMException("AlchemicalPMEForce must define parameters for every particle in the System");
        }
        std::set<std::pair<int, int>> exceptionPairs;
        if (!AlchemicalPME::buildExceptionPairs(force, exceptionPairs)) {
            throw OpenMMException("AlchemicalPMEForce has duplicate exceptions for the same particle pair");
        }

        auto* data = reinterpret_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
        if (data == nullptr || data->contexts.empty()) {
            throw OpenMMException("CUDA AlchemicalPMEForce kernel could not access CudaContext");
        }
        if (data->contexts.size() != 1) {
            throw OpenMMException("CUDA AlchemicalPMEForce currently supports a single CUDA context only");
        }
        cu_ = data->contexts[0];
        cu_->setAsCurrent();
        cu_->addForce(new AlchemicalPMEForceInfo(force));
        cu_->addPostComputation(new AlchemicalPMEDerivativeParameterMirror(context, force));
        cu_->addEnergyParameterDerivative(force.getLambdaParameterName());
        lambdaGlobalParamIndex_ = cu_->registerGlobalParam(force.getLambdaParameterName());
        const auto& derivativeNames = cu_->getEnergyParamDerivNames();
        for (int i = 0; i < static_cast<int>(derivativeNames.size()); ++i) {
            if (derivativeNames[i] == force.getLambdaParameterName()) {
                lambdaDerivativeIndex_ = i;
                break;
            }
        }
        if (lambdaDerivativeIndex_ < 0) {
            throw OpenMMException("AlchemicalPMEForce could not register CUDA lambda energy derivative");
        }

        std::vector<AlchemicalPME::ParticleDeviceParameter> particleParams;
        std::vector<AlchemicalPME::ExceptionDeviceParameter> exceptionParams;
        AlchemicalPME::buildParticleDeviceParameters(force, particleParams);
        AlchemicalPME::buildExceptionDeviceParameters(force, exceptionParams);
        initializeDeviceParameters(particleParams);
        initializeUtilityParameters(particleParams, exceptionParams);
        registerUtilityInteractions(context, force);

        threadBlockSize_ = std::min(256, cu_->getMaxThreadBlockSize());

        std::map<std::string, std::string> defines;
        defines["PADDED_NUM_ATOMS"] = std::to_string(cu_->getPaddedNumAtoms());
        double alpha = 1.0;
        int nx = 1, ny = 1, nz = 1;
        if (force.getNonbondedMethod() == AlchemicalPMEForce::PME) {
            force.getPMEParameters(alpha, nx, ny, nz);
            if (alpha == 0.0 || nx <= 0 || ny <= 0 || nz <= 0) {
                throw OpenMMException("AlchemicalPMEForce PME mode requires explicit non-zero PME parameters");
            }
        }
        defines["M_PI"] = "3.14159265358979323846";
        defines["NUM_ATOMS"] = std::to_string(numAtoms_);
        defines["PME_ORDER"] = std::to_string(force.getPMEInterpolationOrder());
        defines["GRID_SIZE_X"] = std::to_string(nx);
        defines["GRID_SIZE_Y"] = std::to_string(ny);
        defines["GRID_SIZE_Z"] = std::to_string(nz);
        defines["EWALD_ALPHA"] = std::to_string(alpha);
        defines["RECIP_EXP_FACTOR"] = std::to_string(AlchemicalPME::Pi*AlchemicalPME::Pi/(alpha*alpha));
        defines["EPSILON_FACTOR"] = std::to_string(std::sqrt(AlchemicalPME::CoulombConstant));
        const std::string cudaSource =
            std::string(AlchemicalPMECudaSourcePart1) +
            AlchemicalPMECudaSourcePart2 +
            AlchemicalPMECudaSourcePart3 +
            AlchemicalPMECudaSourcePart4;
        module_ = cu_->createModule(cudaSource, defines);
        if (force.getNonbondedMethod() == AlchemicalPMEForce::PME) {
            initializePme(force);
        }
    }

    double execute(ContextImpl& context, bool includeForces, bool includeEnergy) override {
        cu_->setAsCurrent();
        const double lambda = getContextParameter(context, owner_->getLambdaParameterName(), 0.0);
        uploadBoxVectorsIfNeeded();
        if (owner_->getNonbondedMethod() == AlchemicalPMEForce::PME) {
            return executePme(context, lambda, includeForces, includeEnergy);
        }
        return 0.0;
    }

private:
    static double getContextParameter(ContextImpl& context, const std::string& name, double defaultValue) {
        const auto& params = context.getParameters();
        auto found = params.find(name);
        return found == params.end() ? defaultValue : found->second;
    }

    static std::string realLiteral(double value) {
        return std::to_string(value);
    }

    double getEwaldAlpha() const {
        double alpha = 0.0;
        int nx = 0, ny = 0, nz = 0;
        owner_->getPMEParameters(alpha, nx, ny, nz);
        return alpha;
    }

    struct ParticleDeviceParameterFloat {
        float charge0;
        float chargeDelta;
        float sigma0;
        float sigmaDelta;
        float epsilon0;
        float epsilonDelta;
    };

    const AlchemicalPMEForce* owner_;
    CudaContext* cu_;
    CudaArray particleParams_;
    CudaArray utilityParticleParams1_;
    CudaArray utilityParticleParams2_;
    CudaArray utilityExceptionAtoms_;
    CudaArray utilityExceptionParams1_;
    CudaArray utilityExceptionParams2_;
    CudaArray utilityBoxVectors_;
    CudaArray pmeGrid1_;
    CudaArray pmeGrid2_;
    CudaArray pmeBsplineModuliX_;
    CudaArray pmeBsplineModuliY_;
    CudaArray pmeBsplineModuliZ_;
    CudaArray pmeAtomGridIndex_;
    CudaArray pmeCharges_;
    CudaArray pmeDCharges_;
    CudaArray pmeChargeDerivatives_;
    CudaArray pmeEnergyBuffer_;
    CudaArray pmeScalarBuffer_;
    FFT3D pmeFft_;
    CUmodule module_;
    CUfunction pmeGridIndexKernel_;
    CUfunction pmeSpreadChargeKernel_;
    CUfunction pmeConvolutionKernel_;
    CUfunction pmeEnergyKernel_;
    CUfunction pmeInterpolateKernel_;
    CUfunction pmeChargeAndScalarKernel_;
    CUfunction pmeDerivativeReductionKernel_;
    CUfunction pmeFinishKernel_;
    int numAtoms_;
    int numExceptions_;
    int lambdaDerivativeIndex_;
    int lambdaGlobalParamIndex_;
    int threadBlockSize_;
    int pmeGridSizeX_;
    int pmeGridSizeY_;
    int pmeGridSizeZ_;
    int pmeGridElements_;
    int pmeWorkSize_;
    bool cachedBoxVectorsValid_ = false;
    Vec3 cachedBoxVectors_[3];

    double executePme(ContextImpl& context, double lambda, bool includeForces, bool includeEnergy) {
        executeGpuPmeReciprocal(lambda, includeForces, includeEnergy);
        return 0.0;
    }

    void initializePme(const AlchemicalPMEForce& force) {
        double alpha;
        force.getPMEParameters(alpha, pmeGridSizeX_, pmeGridSizeY_, pmeGridSizeZ_);
        pmeGridElements_ = pmeGridSizeX_*pmeGridSizeY_*pmeGridSizeZ_;
        pmeWorkSize_ = ((pmeGridElements_ + threadBlockSize_ - 1)/threadBlockSize_)*threadBlockSize_;
        const int realSize = cu_->getUseDoublePrecision() ? sizeof(double) : sizeof(float);
        const int mixedSize = (cu_->getUseDoublePrecision() || cu_->getUseMixedPrecision()) ? sizeof(double) : sizeof(float);
        pmeGrid1_.initialize(*cu_, pmeGridElements_, 2*realSize, "alchemicalPMEPmeGrid1");
        pmeGrid2_.initialize(*cu_, pmeGridElements_, 2*realSize, "alchemicalPMEPmeGrid2");
        pmeBsplineModuliX_.initialize(*cu_, pmeGridSizeX_, realSize, "alchemicalPMEBsplineX");
        pmeBsplineModuliY_.initialize(*cu_, pmeGridSizeY_, realSize, "alchemicalPMEBsplineY");
        pmeBsplineModuliZ_.initialize(*cu_, pmeGridSizeZ_, realSize, "alchemicalPMEBsplineZ");
        pmeAtomGridIndex_.initialize(*cu_, std::max(1, numAtoms_), sizeof(int2), "alchemicalPMEAtomGridIndex");
        pmeCharges_.initialize(*cu_, std::max(1, numAtoms_), realSize, "alchemicalPMECharges");
        pmeDCharges_.initialize(*cu_, std::max(1, numAtoms_), realSize, "alchemicalPMEDCharges");
        pmeChargeDerivatives_.initialize(*cu_, std::max(1, numAtoms_), realSize, "alchemicalPMEChargeDerivatives");
        pmeEnergyBuffer_.initialize(*cu_, pmeWorkSize_, mixedSize, "alchemicalPMEPmeEnergy");
        pmeScalarBuffer_.initialize(*cu_, 3, mixedSize, "alchemicalPMEPmeScalarTerms");
        cu_->addAutoclearBuffer(pmeGrid1_);
        cu_->addAutoclearBuffer(pmeGrid2_);
        cu_->addAutoclearBuffer(pmeChargeDerivatives_);
        cu_->addAutoclearBuffer(pmeEnergyBuffer_);
        pmeFft_ = cu_->createFFT(pmeGridSizeX_, pmeGridSizeY_, pmeGridSizeZ_, true);

        uploadBsplineModuli(pmeGridSizeX_, pmeBsplineModuliX_);
        uploadBsplineModuli(pmeGridSizeY_, pmeBsplineModuliY_);
        uploadBsplineModuli(pmeGridSizeZ_, pmeBsplineModuliZ_);

        pmeGridIndexKernel_ = cu_->getKernel(module_, "findAlchemicalPMEAtomGridIndex");
        pmeSpreadChargeKernel_ = cu_->getKernel(module_, "spreadAlchemicalPMECharge");
        pmeConvolutionKernel_ = cu_->getKernel(module_, "convolveAlchemicalPMEGrid");
        pmeEnergyKernel_ = cu_->getKernel(module_, "evaluateAlchemicalPMEReciprocalEnergy");
        pmeInterpolateKernel_ = cu_->getKernel(module_, "interpolateAlchemicalPMEForceAndDerivative");
        pmeChargeAndScalarKernel_ = cu_->getKernel(module_, "computeAlchemicalPMEChargesAndScalarTerms");
        pmeDerivativeReductionKernel_ = cu_->getKernel(module_, "reduceAlchemicalPMEReciprocalDerivative");
        pmeFinishKernel_ = cu_->getKernel(module_, "finishAlchemicalPMEDeviceBuffers");
    }

    void uploadRealArray(const std::vector<double>& values, CudaArray& array) {
        if (cu_->getUseDoublePrecision()) {
            array.upload(values);
        }
        else {
            std::vector<float> valuesFloat;
            valuesFloat.reserve(values.size());
            for (double value : values) {
                valuesFloat.push_back(static_cast<float>(value));
            }
            array.upload(valuesFloat);
        }
    }

    void initializeDeviceParameters(const std::vector<AlchemicalPME::ParticleDeviceParameter>& particleParams) {
        const int realSize = cu_->getUseDoublePrecision() ? sizeof(double) : sizeof(float);
        particleParams_.initialize(*cu_, std::max(1, numAtoms_), 6*realSize,
                                   "alchemicalPMEParticleParams");
        if (cu_->getUseDoublePrecision()) {
            if (!particleParams.empty()) {
                particleParams_.upload(particleParams);
            }
        }
        else {
            std::vector<ParticleDeviceParameterFloat> particleBuffer(particleParams.size());
            for (size_t i = 0; i < particleParams.size(); ++i) {
                const auto& p = particleParams[i];
                particleBuffer[i] = {
                    static_cast<float>(p.charge0),
                    static_cast<float>(p.chargeDelta),
                    static_cast<float>(p.sigma0),
                    static_cast<float>(p.sigmaDelta),
                    static_cast<float>(p.epsilon0),
                    static_cast<float>(p.epsilonDelta)
                };
            }
            if (!particleBuffer.empty()) {
                particleParams_.upload(particleBuffer);
            }
        }
    }

    void uploadBsplineModuli(int ndata, CudaArray& array) {
        const int order = owner_->getPMEInterpolationOrder();
        std::vector<double> data(order);
        std::vector<double> ddata(order);
        std::vector<double> bsplinesData(ndata, 0.0);
        data[order-1] = 0.0;
        data[1] = 0.0;
        data[0] = 1.0;
        for (int k = 3; k < order; k++) {
            double div = 1.0/(k-1.0);
            data[k-1] = 0.0;
            for (int l = 1; l < (k-1); l++) {
                data[k-l-1] = div*(l*data[k-l-2] + (k-l)*data[k-l-1]);
            }
            data[0] = div*data[0];
        }
        ddata[0] = -data[0];
        for (int k = 1; k < order; k++) {
            ddata[k] = data[k-1] - data[k];
        }
        double div = 1.0/(order-1);
        data[order-1] = 0.0;
        for (int l = 1; l < (order-1); l++) {
            data[order-l-1] = div*(l*data[order-l-2] + (order-l)*data[order-l-1]);
        }
        data[0] = div*data[0];
        for (int i = 1; i <= order && i < ndata; i++) {
            bsplinesData[i] = data[i-1];
        }
        std::vector<double> moduli(ndata);
        for (int i = 0; i < ndata; i++) {
            double sc = 0.0;
            double ss = 0.0;
            for (int j = 0; j < ndata; j++) {
                double arg = (2.0*AlchemicalPME::Pi*i*j)/ndata;
                sc += bsplinesData[j]*std::cos(arg);
                ss += bsplinesData[j]*std::sin(arg);
            }
            moduli[i] = sc*sc + ss*ss;
        }
        for (int i = 0; i < ndata; i++) {
            if (moduli[i] < 1.0e-7) {
                moduli[i] = 0.5*(moduli[(i-1+ndata)%ndata] + moduli[(i+1)%ndata]);
            }
        }
        uploadRealArray(moduli, array);
    }

    void executeGpuPmeReciprocal(double lambda, bool includeForces, bool includeEnergy) {
        cu_->clearBuffer(pmeGrid1_);
        cu_->clearBuffer(pmeGrid2_);
        cu_->clearBuffer(pmeChargeDerivatives_);
        cu_->clearBuffer(pmeEnergyBuffer_);
        cu_->clearBuffer(pmeScalarBuffer_);
        void* boxSize = cu_->getPeriodicBoxSizePointer();
        void* invBoxSize = cu_->getInvPeriodicBoxSizePointer();
        void* boxX = cu_->getPeriodicBoxVecXPointer();
        void* boxY = cu_->getPeriodicBoxVecYPointer();
        void* boxZ = cu_->getPeriodicBoxVecZPointer();
        Vec3 boxVectors[3];
        cu_->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        Real4Args recipArgs;
        setPmeReciprocalArgs(recipArgs, boxVectors);
        double alpha = 0.0;
        int nx = 0, ny = 0, nz = 0;
        owner_->getPMEParameters(alpha, nx, ny, nz);
        double electrostaticLambdaStart = owner_->getElectrostaticLambdaStart();
        float lambdaFloat = static_cast<float>(lambda);
        float electrostaticLambdaStartFloat = static_cast<float>(electrostaticLambdaStart);
        float alphaFloat = static_cast<float>(alpha);
        void* lambdaArg = cu_->getUseDoublePrecision() ? static_cast<void*>(&lambda) : static_cast<void*>(&lambdaFloat);
        void* electrostaticLambdaStartArg = cu_->getUseDoublePrecision() ? static_cast<void*>(&electrostaticLambdaStart) : static_cast<void*>(&electrostaticLambdaStartFloat);
        void* alphaArg = cu_->getUseDoublePrecision() ? static_cast<void*>(&alpha) : static_cast<void*>(&alphaFloat);
        const int realSize = cu_->getUseDoublePrecision() ? sizeof(double) : sizeof(float);
        const int mixedSize = (cu_->getUseDoublePrecision() || cu_->getUseMixedPrecision()) ? sizeof(double) : sizeof(float);
        void* argsCharges[] = {
            &particleParams_.getDevicePointer(),
            lambdaArg,
            electrostaticLambdaStartArg,
            alphaArg,
            boxX,
            boxY,
            boxZ,
            &pmeCharges_.getDevicePointer(),
            &pmeDCharges_.getDevicePointer(),
            &pmeScalarBuffer_.getDevicePointer()
        };
        cu_->executeKernel(pmeChargeAndScalarKernel_, argsCharges, threadBlockSize_, threadBlockSize_,
                           4*threadBlockSize_*realSize);
        void* argsIndex[] = {
            &cu_->getPosq().getDevicePointer(),
            &pmeAtomGridIndex_.getDevicePointer(),
            boxSize, invBoxSize, boxX, boxY, boxZ,
            recipArgs.xPointer(), recipArgs.yPointer(), recipArgs.zPointer()
        };
        cu_->executeKernel(pmeGridIndexKernel_, argsIndex, numAtoms_, threadBlockSize_);
        void* argsSpread[] = {
            &cu_->getPosq().getDevicePointer(),
            &pmeGrid1_.getDevicePointer(),
            boxSize, invBoxSize, boxX, boxY, boxZ,
            recipArgs.xPointer(), recipArgs.yPointer(), recipArgs.zPointer(),
            &pmeAtomGridIndex_.getDevicePointer(),
            &cu_->getAtomIndexArray().getDevicePointer(),
            &pmeCharges_.getDevicePointer()
        };
        cu_->executeKernel(pmeSpreadChargeKernel_, argsSpread,
                           numAtoms_*owner_->getPMEInterpolationOrder(), threadBlockSize_);
        pmeFft_->execFFT(pmeGrid1_, pmeGrid2_, true);

        if (includeEnergy) {
            void* argsEnergy[] = {
                &pmeGrid2_.getDevicePointer(),
                &pmeBsplineModuliX_.getDevicePointer(),
                &pmeBsplineModuliY_.getDevicePointer(),
                &pmeBsplineModuliZ_.getDevicePointer(),
                recipArgs.xPointer(), recipArgs.yPointer(), recipArgs.zPointer(),
                &pmeEnergyBuffer_.getDevicePointer()
            };
            cu_->executeKernel(pmeEnergyKernel_, argsEnergy, pmeWorkSize_, threadBlockSize_);
        }

        void* argsConvolution[] = {
            &pmeGrid2_.getDevicePointer(),
            &pmeBsplineModuliX_.getDevicePointer(),
            &pmeBsplineModuliY_.getDevicePointer(),
            &pmeBsplineModuliZ_.getDevicePointer(),
            recipArgs.xPointer(), recipArgs.yPointer(), recipArgs.zPointer()
        };
        cu_->executeKernel(pmeConvolutionKernel_, argsConvolution, pmeGridSizeX_*pmeGridSizeY_*pmeGridSizeZ_, threadBlockSize_);
        pmeFft_->execFFT(pmeGrid2_, pmeGrid1_, false);
        int includeForcesInt = includeForces ? 1 : 0;
        void* argsInterpolate[] = {
            &cu_->getPosq().getDevicePointer(),
            &cu_->getForce().getDevicePointer(),
            &pmeGrid1_.getDevicePointer(),
            boxSize, invBoxSize, boxX, boxY, boxZ,
            recipArgs.xPointer(), recipArgs.yPointer(), recipArgs.zPointer(),
            &pmeAtomGridIndex_.getDevicePointer(),
            &cu_->getAtomIndexArray().getDevicePointer(),
            &pmeCharges_.getDevicePointer(),
            &includeForcesInt,
            &pmeChargeDerivatives_.getDevicePointer()
        };
        cu_->executeKernel(pmeInterpolateKernel_, argsInterpolate, numAtoms_, threadBlockSize_);
        void* argsDeriv[] = {
            &pmeChargeDerivatives_.getDevicePointer(),
            &pmeDCharges_.getDevicePointer(),
            &pmeScalarBuffer_.getDevicePointer()
        };
        cu_->executeKernel(pmeDerivativeReductionKernel_, argsDeriv, threadBlockSize_, threadBlockSize_,
                           threadBlockSize_*mixedSize);
        int pmeWorkSize = pmeWorkSize_;
        int energyBufferSize = cu_->getEnergyBuffer().getSize();
        int numDerivatives = static_cast<int>(cu_->getEnergyParamDerivNames().size());
        int includeEnergyInt = includeEnergy ? 1 : 0;
        void* argsFinish[] = {
            &pmeEnergyBuffer_.getDevicePointer(),
            &pmeScalarBuffer_.getDevicePointer(),
            &cu_->getEnergyBuffer().getDevicePointer(),
            &cu_->getEnergyParamDerivBuffer().getDevicePointer(),
            &pmeWorkSize,
            &energyBufferSize,
            &numDerivatives,
            &lambdaDerivativeIndex_,
            &includeEnergyInt
        };
        cu_->executeKernel(pmeFinishKernel_, argsFinish, energyBufferSize, threadBlockSize_);
    }

    struct Real4Args {
        double4 xDouble;
        double4 yDouble;
        double4 zDouble;
        float4 xFloat;
        float4 yFloat;
        float4 zFloat;
        bool useDouble = true;
        void* xPointer() {
            return useDouble ? static_cast<void*>(&xDouble) : static_cast<void*>(&xFloat);
        }
        void* yPointer() {
            return useDouble ? static_cast<void*>(&yDouble) : static_cast<void*>(&yFloat);
        }
        void* zPointer() {
            return useDouble ? static_cast<void*>(&zDouble) : static_cast<void*>(&zFloat);
        }
    };

    void setReal4Args(Real4Args& values, const Vec3 recip[3]) {
        values.useDouble = cu_->getUseDoublePrecision();
        values.xDouble = make_double4(recip[0][0], recip[0][1], recip[0][2], 0.0);
        values.yDouble = make_double4(recip[1][0], recip[1][1], recip[1][2], 0.0);
        values.zDouble = make_double4(recip[2][0], recip[2][1], recip[2][2], 0.0);
        values.xFloat = make_float4(static_cast<float>(recip[0][0]), static_cast<float>(recip[0][1]), static_cast<float>(recip[0][2]), 0.0f);
        values.yFloat = make_float4(static_cast<float>(recip[1][0]), static_cast<float>(recip[1][1]), static_cast<float>(recip[1][2]), 0.0f);
        values.zFloat = make_float4(static_cast<float>(recip[2][0]), static_cast<float>(recip[2][1]), static_cast<float>(recip[2][2]), 0.0f);
    }

    void setPmeReciprocalArgs(Real4Args& values, const Vec3 boxVectors[3]) {
        const double determinant = boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2];
        if (std::abs(determinant) < std::numeric_limits<double>::epsilon()) {
            throw OpenMMException("AlchemicalPMEForce requires a non-singular periodic box");
        }
        const double scale = 1.0/determinant;
        Vec3 pmeRecip[3];
        pmeRecip[0] = Vec3(boxVectors[1][1]*boxVectors[2][2]*scale, 0.0, 0.0);
        pmeRecip[1] = Vec3(-boxVectors[1][0]*boxVectors[2][2]*scale,
                           boxVectors[0][0]*boxVectors[2][2]*scale, 0.0);
        pmeRecip[2] = Vec3((boxVectors[1][0]*boxVectors[2][1] - boxVectors[1][1]*boxVectors[2][0])*scale,
                           -boxVectors[0][0]*boxVectors[2][1]*scale,
                           boxVectors[0][0]*boxVectors[1][1]*scale);
        setReal4Args(values, pmeRecip);
    }

    void uploadReal4Array(const std::vector<std::array<double, 4>>& values, CudaArray& array) {
        if (cu_->getUseDoublePrecision()) {
            std::vector<double4> buffer(values.size());
            for (size_t i = 0; i < values.size(); ++i) {
                buffer[i] = make_double4(values[i][0], values[i][1], values[i][2], values[i][3]);
            }
            array.upload(buffer);
        }
        else {
            std::vector<float4> buffer(values.size());
            for (size_t i = 0; i < values.size(); ++i) {
                buffer[i] = make_float4(static_cast<float>(values[i][0]),
                                        static_cast<float>(values[i][1]),
                                        static_cast<float>(values[i][2]),
                                        static_cast<float>(values[i][3]));
            }
            array.upload(buffer);
        }
    }

    void uploadBoxVectorsIfNeeded() {
        Vec3 boxVectors[3];
        cu_->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        bool changed = !cachedBoxVectorsValid_;
        for (int i = 0; i < 3 && !changed; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (cachedBoxVectors_[i][j] != boxVectors[i][j]) {
                    changed = true;
                    break;
                }
            }
        }
        if (!changed) {
            return;
        }
        std::vector<std::array<double, 4>> boxData(4);
        boxData[0] = {boxVectors[0][0], boxVectors[0][1], boxVectors[0][2], 0.0};
        boxData[1] = {boxVectors[1][0], boxVectors[1][1], boxVectors[1][2], 0.0};
        boxData[2] = {boxVectors[2][0], boxVectors[2][1], boxVectors[2][2], 0.0};
        boxData[3] = {
            boxVectors[0][0] == 0.0 ? 0.0 : 1.0/boxVectors[0][0],
            boxVectors[1][1] == 0.0 ? 0.0 : 1.0/boxVectors[1][1],
            boxVectors[2][2] == 0.0 ? 0.0 : 1.0/boxVectors[2][2],
            0.0
        };
        uploadReal4Array(boxData, utilityBoxVectors_);
        for (int i = 0; i < 3; ++i) {
            cachedBoxVectors_[i] = boxVectors[i];
        }
        cachedBoxVectorsValid_ = true;
    }

    void initializeUtilityParameters(const std::vector<AlchemicalPME::ParticleDeviceParameter>& particleParams,
                                     const std::vector<AlchemicalPME::ExceptionDeviceParameter>& exceptionParams) {
        const int realSize = cu_->getUseDoublePrecision() ? sizeof(double) : sizeof(float);
        utilityParticleParams1_.initialize(*cu_, std::max(1, numAtoms_), 4*realSize,
                                           "alchemicalPMEUtilityParticleParams1");
        utilityParticleParams2_.initialize(*cu_, std::max(1, numAtoms_), 4*realSize,
                                           "alchemicalPMEUtilityParticleParams2");
        utilityExceptionAtoms_.initialize<int2>(*cu_, std::max(1, numExceptions_),
                                                "alchemicalPMEUtilityExceptionAtoms");
        utilityExceptionParams1_.initialize(*cu_, std::max(1, numExceptions_), 4*realSize,
                                            "alchemicalPMEUtilityExceptionParams1");
        utilityExceptionParams2_.initialize(*cu_, std::max(1, numExceptions_), 4*realSize,
                                            "alchemicalPMEUtilityExceptionParams2");
        utilityBoxVectors_.initialize(*cu_, 4, 4*realSize, "alchemicalPMEUtilityBoxVectors");

        std::vector<std::array<double, 4>> particle1(std::max(1, numAtoms_));
        std::vector<std::array<double, 4>> particle2(std::max(1, numAtoms_));
        for (int i = 0; i < numAtoms_; ++i) {
            const auto& p = particleParams[i];
            particle1[i] = {p.charge0, p.chargeDelta, p.sigma0, p.sigmaDelta};
            particle2[i] = {p.epsilon0, p.epsilonDelta, 0.0, 0.0};
        }
        uploadReal4Array(particle1, utilityParticleParams1_);
        uploadReal4Array(particle2, utilityParticleParams2_);

        std::vector<int2> exceptionAtoms(std::max(1, numExceptions_));
        std::vector<std::array<double, 4>> exception1(std::max(1, numExceptions_));
        std::vector<std::array<double, 4>> exception2(std::max(1, numExceptions_));
        for (int i = 0; i < numExceptions_; ++i) {
            const auto& e = exceptionParams[i];
            exceptionAtoms[i] = make_int2(e.particle1, e.particle2);
            exception1[i] = {e.chargeProd0, e.chargeProdDelta, e.sigma0, e.sigmaDelta};
            exception2[i] = {e.epsilon0, e.epsilonDelta, 0.0, 0.0};
        }
        utilityExceptionAtoms_.upload(exceptionAtoms);
        uploadReal4Array(exception1, utilityExceptionParams1_);
        uploadReal4Array(exception2, utilityExceptionParams2_);
        uploadBoxVectorsIfNeeded();
    }

    std::string utilityCommonSource(const std::string& derivVariable,
                                    const std::string& lambdaName,
                                    const std::string& particleParams1Name,
                                    const std::string& particleParams2Name,
                                    const std::string& exceptionParams1Name,
                                    const std::string& exceptionParams2Name,
                                    bool pmeException) const {
        const double alpha = getEwaldAlpha();
        const double twoOverSqrtPi = 2.0/AlchemicalPME::SqrtPi;
        std::string source;
        source += "const real alchLambda = "+lambdaName+"[0];\n";
        source += "const real alchElecStart = (real) "+realLiteral(owner_->getElectrostaticLambdaStart())+";\n";
        source += "real alchClamp = (alchLambda < 0 ? 0 : (alchLambda > 1 ? 1 : alchLambda));\n";
        source += "real alchActive = (alchElecStart <= 0 ? alchLambda : (alchLambda-alchElecStart)/(1-alchElecStart));\n";
        source += "alchActive = (alchActive < 0 ? 0 : (alchActive > 1 ? 1 : alchActive));\n";
        source += "real alchDActive = (alchElecStart <= 0 ? 1 : ((alchLambda <= alchElecStart || alchLambda >= 1) ? 0 : RECIP(1-alchElecStart)));\n";
        source += "real tempAlchEnergy = 0;\n";
        source += "real tempAlchDeriv = 0;\n";
        source += "real tempAlchForce = 0;\n";
        source += "real sigma = 0;\n";
        source += "real dSigma = 0;\n";
        source += "real epsilon = 0;\n";
        source += "real dEpsilon = 0;\n";
        source += "real chargeProd = 0;\n";
        source += "real dChargeProd = 0;\n";
        source += "real invR2 = invR*invR;\n";
        source += "real r = RECIP(invR);\n";
        source += "real sr = 0;\n";
        source += "real sr2 = 0;\n";
        source += "real sr6 = 0;\n";
        source += "real sr12 = 0;\n";
        if (pmeException) {
            source += "const real alchAlpha = (real) "+realLiteral(alpha)+";\n";
            source += "const real alchTwoOverSqrtPi = (real) "+realLiteral(twoOverSqrtPi)+";\n";
            source += "(void) alchAlpha;\n";
            source += "(void) alchTwoOverSqrtPi;\n";
        }
        source += "(void) "+particleParams1Name+";\n";
        source += "(void) "+particleParams2Name+";\n";
        source += "(void) "+exceptionParams1Name+";\n";
        source += "(void) "+exceptionParams2Name+";\n";
        source += "(void) alchClamp;\n";
        source += "(void) alchActive;\n";
        source += "(void) alchDActive;\n";
        source += "(void) tempAlchEnergy;\n";
        source += "(void) tempAlchDeriv;\n";
        source += "(void) tempAlchForce;\n";
        source += "(void) sigma;\n";
        source += "(void) dSigma;\n";
        source += "(void) epsilon;\n";
        source += "(void) dEpsilon;\n";
        source += "(void) chargeProd;\n";
        source += "(void) dChargeProd;\n";
        source += "(void) sr;\n";
        source += "(void) sr2;\n";
        source += "(void) sr6;\n";
        source += "(void) sr12;\n";
        source += "(void) "+derivVariable+";\n";
        return source;
    }

    std::string electrostaticInterpolationSource(const std::string& valueName,
                                                 const std::string& derivName,
                                                 const std::string& value0,
                                                 const std::string& valueDelta,
                                                 const std::string& lambda,
                                                 const std::string& lambdaStart) const {
        std::string source;
        source += "{\n";
        source += "real alchValue0 = "+value0+";\n";
        source += "real alchValueDelta = "+valueDelta+";\n";
        source += "real alchValue1 = alchValue0 + alchValueDelta;\n";
        source += "if ("+lambdaStart+" <= 0) {\n";
        source += "    "+valueName+" = alchValue0 + "+lambda+"*alchValueDelta;\n";
        source += "    "+derivName+" = alchValueDelta;\n";
        source += "}\n";
        source += "else {\n";
        source += "    bool alchGrow = fabs(alchValue0) <= (real) 1.0e-12 && fabs(alchValue1) > (real) 1.0e-12;\n";
        source += "    bool alchVanish = fabs(alchValue1) <= (real) 1.0e-12 && fabs(alchValue0) > (real) 1.0e-12;\n";
        source += "    if (alchGrow) {\n";
        source += "        real alchDenom = 1-"+lambdaStart+";\n";
        source += "        if ("+lambda+" <= "+lambdaStart+") { "+valueName+" = 0; "+derivName+" = 0; }\n";
        source += "        else {\n";
        source += "            real alchSchedule = ("+lambda+"-"+lambdaStart+")/alchDenom;\n";
        source += "            alchSchedule = alchSchedule < 0 ? 0 : (alchSchedule > 1 ? 1 : alchSchedule);\n";
        source += "            "+valueName+" = alchValue1*alchSchedule;\n";
        source += "            "+derivName+" = alchSchedule < 1 ? alchValue1/alchDenom : 0;\n";
        source += "        }\n";
        source += "    }\n";
        source += "    else if (alchVanish) {\n";
        source += "        real alchEnd = 1-"+lambdaStart+";\n";
        source += "        if ("+lambda+" >= alchEnd) { "+valueName+" = 0; "+derivName+" = 0; }\n";
        source += "        else {\n";
        source += "            real alchSchedule = 1-"+lambda+"/alchEnd;\n";
        source += "            alchSchedule = alchSchedule < 0 ? 0 : (alchSchedule > 1 ? 1 : alchSchedule);\n";
        source += "            "+valueName+" = alchValue0*alchSchedule;\n";
        source += "            "+derivName+" = alchSchedule > 0 ? -alchValue0/alchEnd : 0;\n";
        source += "        }\n";
        source += "    }\n";
        source += "    else {\n";
        source += "        "+valueName+" = alchValue0 + "+lambda+"*alchValueDelta;\n";
        source += "        "+derivName+" = alchValueDelta;\n";
        source += "    }\n";
        source += "}\n";
        source += "}\n";
        return source;
    }

    std::string buildDirectInteractionSource(const std::string& derivVariable) const {
        const bool pme = owner_->getNonbondedMethod() == AlchemicalPMEForce::PME;
        const bool usesCutoff = owner_->getNonbondedMethod() != AlchemicalPMEForce::NoCutoff;
        const bool reactionField = owner_->getNonbondedMethod() == AlchemicalPMEForce::CutoffPeriodic;
        const double cutoff2 = owner_->getCutoffDistance()*owner_->getCutoffDistance();
        const double alpha = getEwaldAlpha();
        const double twoOverSqrtPi = 2.0/AlchemicalPME::SqrtPi;
        double reactionFieldK = 0.0;
        double reactionFieldC = 0.0;
        if (reactionField) {
            const double cutoff = owner_->getCutoffDistance();
            const double dielectric = owner_->getReactionFieldDielectric();
            reactionFieldK = (dielectric - 1.0)/(2.0*dielectric + 1.0)/(cutoff*cutoff*cutoff);
            reactionFieldC = (3.0*dielectric)/(2.0*dielectric + 1.0)/cutoff;
        }
        std::string source;
        source += usesCutoff ? "if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 > 0 && r2 < (real) "+realLiteral(cutoff2)+") {\n" : "if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && r2 > 0) {\n";
        source += "const real alchLambda = alchemicalPMEUtilityGlobals["+std::to_string(lambdaGlobalParamIndex_)+"];\n";
        source += "const real alchElecStart = (real) "+realLiteral(owner_->getElectrostaticLambdaStart())+";\n";
        source += "const int alchAtom1 = alchemicalPMEAtomIndex[atom1];\n";
        source += "const int alchAtom2 = alchemicalPMEAtomIndex[atom2];\n";
        source += "unsigned int alchIncludeInteraction = !isExcluded;\n";
        source += "if (alchIncludeInteraction) {\n";
        source += "const real4 alchemicalPMEUtilityParticleParams11 = alchemicalPMEUtilityParticleParams1[alchAtom1];\n";
        source += "const real4 alchemicalPMEUtilityParticleParams12 = alchemicalPMEUtilityParticleParams1[alchAtom2];\n";
        source += "const real4 alchemicalPMEUtilityParticleParams21 = alchemicalPMEUtilityParticleParams2[alchAtom1];\n";
        source += "const real4 alchemicalPMEUtilityParticleParams22 = alchemicalPMEUtilityParticleParams2[alchAtom2];\n";
        source += "real alchActive = (alchElecStart <= 0 ? alchLambda : (alchLambda-alchElecStart)/(1-alchElecStart));\n";
        source += "alchActive = (alchActive < 0 ? 0 : (alchActive > 1 ? 1 : alchActive));\n";
        source += "real alchDActive = (alchElecStart <= 0 ? 1 : ((alchLambda <= alchElecStart || alchLambda >= 1) ? 0 : RECIP(1-alchElecStart)));\n";
        source += "real qi = 0;\n";
        source += "real dqi = 0;\n";
        source += electrostaticInterpolationSource("qi", "dqi", "alchemicalPMEUtilityParticleParams11.x",
                                                   "alchemicalPMEUtilityParticleParams11.y",
                                                   "alchLambda", "alchElecStart");
        source += "real qj = 0;\n";
        source += "real dqj = 0;\n";
        source += electrostaticInterpolationSource("qj", "dqj", "alchemicalPMEUtilityParticleParams12.x",
                                                   "alchemicalPMEUtilityParticleParams12.y",
                                                   "alchLambda", "alchElecStart");
        source += "real sigmaI = alchemicalPMEUtilityParticleParams11.z + alchLambda*alchemicalPMEUtilityParticleParams11.w;\n";
        source += "real sigmaJ = alchemicalPMEUtilityParticleParams12.z + alchLambda*alchemicalPMEUtilityParticleParams12.w;\n";
        source += "real dSigma = 0.5f*(alchemicalPMEUtilityParticleParams11.w + alchemicalPMEUtilityParticleParams12.w);\n";
        source += "real epsilonI = alchemicalPMEUtilityParticleParams21.x + alchLambda*alchemicalPMEUtilityParticleParams21.y;\n";
        source += "real epsilonJ = alchemicalPMEUtilityParticleParams22.x + alchLambda*alchemicalPMEUtilityParticleParams22.y;\n";
        source += "real dEpsilonI = alchemicalPMEUtilityParticleParams21.y;\n";
        source += "real dEpsilonJ = alchemicalPMEUtilityParticleParams22.y;\n";
        source += "real chargeProd = qi*qj;\n";
        source += "real dChargeProd = dqi*qj + qi*dqj;\n";
        source += "real invR2 = invR*invR;\n";
        source += "real r = RECIP(invR);\n";
        if (pme) {
            source += "real alphaR = (real) "+realLiteral(alpha)+"*r;\n";
            source += "real expAlphaRSqr = EXP(-alphaR*alphaR);\n";
            source += "real erfcAlphaR = ERFC(alphaR);\n";
            source += "real prefactor = (real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*chargeProd*invR;\n";
            source += "tempEnergy += prefactor*erfcAlphaR;\n";
            source += derivVariable+" += interactionScale*(real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*dChargeProd*invR*erfcAlphaR;\n";
            source += "dEdR += prefactor*(erfcAlphaR + alphaR*expAlphaRSqr*(real) "+realLiteral(twoOverSqrtPi)+")*invR2;\n";
        }
        else {
            source += "real prefactor = (real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*chargeProd;\n";
            if (reactionField) {
                source += "real reactionFieldK = (real) "+realLiteral(reactionFieldK)+";\n";
                source += "real reactionFieldC = (real) "+realLiteral(reactionFieldC)+";\n";
                source += "real reactionFieldEnergy = invR + reactionFieldK*r2 - reactionFieldC;\n";
                source += "tempEnergy += prefactor*reactionFieldEnergy;\n";
                source += derivVariable+" += interactionScale*(real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*dChargeProd*reactionFieldEnergy;\n";
                source += "dEdR += prefactor*(invR2*invR - 2*reactionFieldK);\n";
            }
            else {
                source += "tempEnergy += prefactor*invR;\n";
                source += derivVariable+" += interactionScale*(real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*dChargeProd*invR;\n";
                source += "dEdR += prefactor*invR2*invR;\n";
            }
        }
        source += "real sigma = 0.5f*(sigmaI + sigmaJ);\n";
        source += "real epsilon = SQRT(epsilonI*epsilonJ);\n";
        source += "if (sigma != 0 && epsilon != 0) {\n";
        source += "    real dEpsilon = (epsilonI != 0 && epsilonJ != 0 ? 0.5f*epsilon*(dEpsilonI/epsilonI + dEpsilonJ/epsilonJ) : 0);\n";
        source += "    real sr = sigma*invR;\n";
        source += "    real sr2 = sr*sr;\n";
        source += "    real sr6 = sr2*sr2*sr2;\n";
        source += "    real sr12 = sr6*sr6;\n";
        source += "    tempEnergy += 4*epsilon*(sr12-sr6);\n";
        source += derivVariable+" += interactionScale*(4*dEpsilon*(sr12-sr6) + 4*epsilon*(12*sr12-6*sr6)*(dSigma/sigma));\n";
        source += "    dEdR += 24*epsilon*(2*sr12-sr6)*invR2;\n";
        source += "}\n";
        source += "}\n";
        source += "}\n";
        return source;
    }

    std::string buildExceptionInteractionSource(const std::string& derivVariable,
                                                const std::string& lambdaArg,
                                                const std::string& exceptionAtomsArg,
                                                const std::string& exceptionParams1Arg,
                                                const std::string& exceptionParams2Arg,
                                                const std::string& particleParams1Arg,
                                                const std::string& particleParams2Arg,
                                                int exceptionStartIndex) const {
        const bool pme = owner_->getNonbondedMethod() == AlchemicalPMEForce::PME;
        const bool periodic = owner_->getNonbondedMethod() != AlchemicalPMEForce::NoCutoff;
        std::string source;
        source += "real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);\n";
        if (periodic) {
            source += "APPLY_PERIODIC_TO_DELTA(delta)\n";
        }
        source += "real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;\n";
        source += "real invR = RSQRT(r2);\n";
        source += "real invR2 = invR*invR;\n";
        source += "real r = RECIP(invR);\n";
        source += "real alchLambda = "+lambdaArg+"["+std::to_string(lambdaGlobalParamIndex_)+"];\n";
        source += "real alchElecStart = (real) "+realLiteral(owner_->getElectrostaticLambdaStart())+";\n";
        source += "real alchActive = (alchElecStart <= 0 ? alchLambda : (alchLambda-alchElecStart)/(1-alchElecStart));\n";
        source += "alchActive = (alchActive < 0 ? 0 : (alchActive > 1 ? 1 : alchActive));\n";
        source += "real alchDActive = (alchElecStart <= 0 ? 1 : ((alchLambda <= alchElecStart || alchLambda >= 1) ? 0 : RECIP(1-alchElecStart)));\n";
        source += "int exceptionIndex = index + "+std::to_string(exceptionStartIndex)+";\n";
        source += "int2 exceptionAtoms = "+exceptionAtomsArg+"[index];\n";
        source += "(void) exceptionAtoms;\n";
        source += "real4 exceptionParams1 = "+exceptionParams1Arg+"[exceptionIndex];\n";
        source += "real4 exceptionParams2 = "+exceptionParams2Arg+"[exceptionIndex];\n";
        source += "real exceptionChargeProd = 0;\n";
        source += "real exceptionDChargeProd = 0;\n";
        source += electrostaticInterpolationSource("exceptionChargeProd", "exceptionDChargeProd",
                                                   "exceptionParams1.x", "exceptionParams1.y",
                                                   "alchLambda", "alchElecStart");
        source += "real sigma = exceptionParams1.z + alchLambda*exceptionParams1.w;\n";
        source += "real dSigma = exceptionParams1.w;\n";
        source += "real epsilon = exceptionParams2.x + alchLambda*exceptionParams2.y;\n";
        source += "real dEpsilon = exceptionParams2.y;\n";
        source += "real tempEnergy = 0;\n";
        source += "real tempForce = 0;\n";
        source += "tempEnergy += (real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*exceptionChargeProd*invR;\n";
        source += derivVariable+" += (real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*exceptionDChargeProd*invR;\n";
        source += "tempForce += (real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*exceptionChargeProd*invR;\n";
        source += "if (sigma != 0 && epsilon != 0) {\n";
        source += "    real sr = sigma*invR;\n";
        source += "    real sr2 = sr*sr;\n";
        source += "    real sr6 = sr2*sr2*sr2;\n";
        source += "    real sr12 = sr6*sr6;\n";
        source += "    tempEnergy += 4*epsilon*(sr12-sr6);\n";
        source += derivVariable+" += 4*dEpsilon*(sr12-sr6) + 4*epsilon*(12*sr12-6*sr6)*(dSigma/sigma);\n";
        source += "    tempForce += 24*epsilon*(2*sr12-sr6);\n";
        source += "}\n";
        source += "energy += tempEnergy;\n";
        source += "delta *= tempForce*invR2;\n";
        source += "real3 force1 = -delta;\n";
        source += "real3 force2 = delta;\n";
        return source;
    }

    std::string buildPmeExclusionInteractionSource(const std::string& derivVariable,
                                                   const std::string& lambdaArg,
                                                   const std::string& exceptionAtomsArg,
                                                   const std::string& particleParams1Arg) const {
        const double alpha = getEwaldAlpha();
        const double twoOverSqrtPi = 2.0/AlchemicalPME::SqrtPi;
        std::string source;
        source += "int2 exceptionAtoms = "+exceptionAtomsArg+"[index];\n";
        source += "real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);\n";
        source += "APPLY_PERIODIC_TO_DELTA(delta)\n";
        source += "real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;\n";
        source += "real invR = RSQRT(r2);\n";
        source += "real r = RECIP(invR);\n";
        source += "real alchLambda = "+lambdaArg+"["+std::to_string(lambdaGlobalParamIndex_)+"];\n";
        source += "real alchElecStart = (real) "+realLiteral(owner_->getElectrostaticLambdaStart())+";\n";
        source += "real4 particle1Params1 = "+particleParams1Arg+"[exceptionAtoms.x];\n";
        source += "real4 particle2Params1 = "+particleParams1Arg+"[exceptionAtoms.y];\n";
        source += "real qi = 0;\n";
        source += "real dqi = 0;\n";
        source += electrostaticInterpolationSource("qi", "dqi", "particle1Params1.x",
                                                   "particle1Params1.y", "alchLambda", "alchElecStart");
        source += "real qj = 0;\n";
        source += "real dqj = 0;\n";
        source += electrostaticInterpolationSource("qj", "dqj", "particle2Params1.x",
                                                   "particle2Params1.y", "alchLambda", "alchElecStart");
        source += "real normalChargeProd = qi*qj;\n";
        source += "real normalDChargeProd = dqi*qj + qi*dqj;\n";
        source += "real alphaR = (real) "+realLiteral(alpha)+"*r;\n";
        source += "real erfAlphaR = ERF(alphaR);\n";
        source += "real expAlphaRSqr = EXP(-alphaR*alphaR);\n";
        source += "real tempEnergy = -(real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*normalChargeProd*erfAlphaR*invR;\n";
        source += derivVariable+" -= (real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*normalDChargeProd*erfAlphaR*invR;\n";
        source += "real tempForce = -(real) "+realLiteral(AlchemicalPME::CoulombConstant)+"*normalChargeProd*invR*(erfAlphaR - alphaR*expAlphaRSqr*(real) "+realLiteral(twoOverSqrtPi)+");\n";
        source += "energy += tempEnergy;\n";
        source += "delta *= tempForce*invR*invR;\n";
        source += "real3 force1 = -delta;\n";
        source += "real3 force2 = delta;\n";
        return source;
    }

    void registerUtilityInteractions(ContextImpl& context, const AlchemicalPMEForce& force) {
        std::vector<std::vector<int>> exclusions(numAtoms_);
        std::vector<std::vector<int>> exceptionAtoms(numExceptions_, std::vector<int>(2));
        for (int atom = 0; atom < numAtoms_; ++atom) {
            exclusions[atom].push_back(atom);
        }
        for (int i = 0; i < numExceptions_; ++i) {
            int particle1, particle2;
            double chargeProd0, chargeProd1, sigma0, sigma1, epsilon0, epsilon1;
            force.getExceptionParameters(i, particle1, particle2, chargeProd0, chargeProd1, sigma0, sigma1,
                                         epsilon0, epsilon1);
            exclusions[particle1].push_back(particle2);
            exclusions[particle2].push_back(particle1);
            exceptionAtoms[i][0] = particle1;
            exceptionAtoms[i][1] = particle2;
        }

        auto& nonbonded = cu_->getNonbondedUtilities();
        nonbonded.addArgument(ComputeParameterInfo(utilityParticleParams1_, "alchemicalPMEUtilityParticleParams1",
                                                   "real", 4));
        nonbonded.addArgument(ComputeParameterInfo(utilityParticleParams2_, "alchemicalPMEUtilityParticleParams2",
                                                   "real", 4));
        nonbonded.addArgument(ComputeParameterInfo(cu_->getAtomIndexArray(), "alchemicalPMEAtomIndex",
                                                   "int", 1, false));
        nonbonded.addArgument(ComputeParameterInfo(cu_->getGlobalParamValues(), "alchemicalPMEUtilityGlobals",
                                                   "real", 1));
        nonbonded.addArgument(ComputeParameterInfo(utilityBoxVectors_, "alchemicalPMEUtilityBoxVectors", "real", 4,
                                                   false));
        const std::string nonbondedDeriv = nonbonded.addEnergyParameterDerivative(force.getLambdaParameterName());
        const bool usesCutoff = force.getNonbondedMethod() != AlchemicalPMEForce::NoCutoff;
        const bool usesPeriodic = force.usesPeriodicBoundaryConditions();
        const bool registerDirectInteractions = true;
        if (registerDirectInteractions) {
            nonbonded.addInteraction(usesCutoff, usesPeriodic, numExceptions_ > 0, force.getCutoffDistance(),
                                     exclusions, buildDirectInteractionSource(nonbondedDeriv), force.getForceGroup(),
                                     false, false);
        }

        const bool registerExceptionInteractions = true;
        if (registerExceptionInteractions && numExceptions_ > 0) {
            auto& bonded = cu_->getBondedUtilities();
            const std::string exceptionAtomsArg = bonded.addArgument(utilityExceptionAtoms_, "int2");
            const std::string exceptionParams1Arg = bonded.addArgument(utilityExceptionParams1_, "real4");
            const std::string exceptionParams2Arg = bonded.addArgument(utilityExceptionParams2_, "real4");
            const std::string particleParams1Arg = bonded.addArgument(utilityParticleParams1_, "real4");
            const std::string particleParams2Arg = bonded.addArgument(utilityParticleParams2_, "real4");
            const std::string lambdaArg = bonded.addArgument(cu_->getGlobalParamValues(), "real");
            const std::string bondedDeriv = bonded.addEnergyParameterDerivative(force.getLambdaParameterName());
            if (force.getNonbondedMethod() == AlchemicalPMEForce::PME) {
                bonded.addInteraction(exceptionAtoms,
                                      buildPmeExclusionInteractionSource(bondedDeriv, lambdaArg,
                                                                         exceptionAtomsArg, particleParams1Arg),
                                      force.getForceGroup());
            }
            bonded.addInteraction(exceptionAtoms,
                                  buildExceptionInteractionSource(bondedDeriv, lambdaArg, exceptionAtomsArg,
                                                                  exceptionParams1Arg, exceptionParams2Arg,
                                                                  particleParams1Arg, particleParams2Arg,
                                                                  0),
                                  force.getForceGroup());
        }
    }
};

#else

class CudaCalcAlchemicalPMEForceKernel final : public CalcAlchemicalPMEForceKernel {
public:
    CudaCalcAlchemicalPMEForceKernel(std::string name, const Platform& platform)
        : CalcAlchemicalPMEForceKernel(name, platform) {
    }

    void initialize(ContextImpl& context, const AlchemicalPMEForce& force) override {
        throw OpenMMException(
            "CUDA AlchemicalPMEForce kernel was not built. Reconfigure with "
            "-DOPENMM_COLVARS_ENABLE_CUDA_KERNELS=ON and provide OpenMM CUDA "
            "internal headers.");
    }

    double execute(ContextImpl& context, bool includeForces, bool includeEnergy) override {
        throw OpenMMException("CUDA AlchemicalPMEForce kernel was not built");
    }
};

#endif

class CudaAlchemicalPMEKernelFactory final : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const override {
        if (name == CalcAlchemicalPMEForceKernel::Name()) {
            return new CudaCalcAlchemicalPMEForceKernel(name, platform);
        }
        throw OpenMMException("Tried to create an illegal CUDA AlchemicalPMEForce kernel");
    }
};

} // namespace

void OpenMM::registerCudaAlchemicalPMEKernelFactory(Platform& platform)
{
    platform.registerKernelFactory(CalcAlchemicalPMEForceKernel::Name(),
                                   new CudaAlchemicalPMEKernelFactory());
}
