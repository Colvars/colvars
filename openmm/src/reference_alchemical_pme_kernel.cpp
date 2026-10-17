#include "openmm/internal/AlchemicalPMEKernels.h"

#include "ReferencePME.h"
#include "ReferencePlatform.h"
#include "openmm/KernelFactory.h"
#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "openmm/internal/AlchemicalPMEUtilities.h"
#include "openmm/internal/ContextImpl.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <utility>
#include <vector>

using namespace OpenMM;

namespace {

using AlchemicalPME::CoulombConstant;
using AlchemicalPME::ExceptionState;
using AlchemicalPME::ParticleState;
using AlchemicalPME::Pi;
using AlchemicalPME::SqrtPi;
using AlchemicalPME::buildExceptionPairs;
using AlchemicalPME::buildExceptionStates;
using AlchemicalPME::buildParticleStates;
using AlchemicalPME::cross;
using AlchemicalPME::dot;
using AlchemicalPME::orderedPair;

void reciprocalBoxVectors(const Vec3 boxVectors[3], Vec3 recipBoxVectors[3])
{
    const double volume = dot(boxVectors[0], cross(boxVectors[1], boxVectors[2]));
    if (std::abs(volume) < std::numeric_limits<double>::epsilon()) {
        throw OpenMMException("AlchemicalPMEForce requires a non-singular periodic box");
    }
    recipBoxVectors[0] = cross(boxVectors[1], boxVectors[2])*(1.0/volume);
    recipBoxVectors[1] = cross(boxVectors[2], boxVectors[0])*(1.0/volume);
    recipBoxVectors[2] = cross(boxVectors[0], boxVectors[1])*(1.0/volume);
}

Vec3 minimumImage(const Vec3& delta, const Vec3 boxVectors[3])
{
    Vec3 recip[3];
    reciprocalBoxVectors(boxVectors, recip);
    double x = dot(delta, recip[0]);
    double y = dot(delta, recip[1]);
    double z = dot(delta, recip[2]);
    x -= std::floor(x + 0.5);
    y -= std::floor(y + 0.5);
    z -= std::floor(z + 0.5);
    return boxVectors[0]*x + boxVectors[1]*y + boxVectors[2]*z;
}

double getContextParameter(ContextImpl& context, const std::string& name, double defaultValue)
{
    const auto& params = context.getParameters();
    auto found = params.find(name);
    return found == params.end() ? defaultValue : found->second;
}

class ReferenceCalcAlchemicalPMEForceKernel final : public CalcAlchemicalPMEForceKernel {
public:
    ReferenceCalcAlchemicalPMEForceKernel(std::string name, const Platform& platform)
        : CalcAlchemicalPMEForceKernel(name, platform),
          owner_(nullptr) {
    }

    void initialize(ContextImpl& context, const AlchemicalPMEForce& force) override {
        owner_ = &force;
        if (force.getNumParticles() != context.getSystem().getNumParticles()) {
            throw OpenMMException("AlchemicalPMEForce must define parameters for every particle in the System");
        }
        if (!buildExceptionPairs(force, exceptionPairs_)) {
            throw OpenMMException("AlchemicalPMEForce contains duplicate exceptions");
        }
        localForces_.resize(context.getSystem().getNumParticles());
    }

    double execute(ContextImpl& context, bool includeForces, bool includeEnergy) override {
        ReferencePlatform::PlatformData* data =
            reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
        std::fill(localForces_.begin(), localForces_.end(), Vec3());
        double dEdLambda = 0.0;
        const double energy = calculate(context, *data->positions, localForces_,
                                        data->periodicBoxVectors, dEdLambda);
        if (includeForces) {
            for (int i = 0; i < static_cast<int>(localForces_.size()); ++i) {
                (*data->forces)[i] += localForces_[i];
            }
        }
        if (data->energyParameterDerivatives != nullptr) {
            (*data->energyParameterDerivatives)[owner_->getLambdaParameterName()] += dEdLambda;
        }
        auto& parameters = const_cast<std::map<std::string, double>&>(context.getParameters());
        auto derivativeParameter = parameters.find(owner_->getLambdaDerivativeParameterName());
        if (derivativeParameter != parameters.end()) {
            derivativeParameter->second = dEdLambda;
        }
        return includeEnergy ? energy : 0.0;
    }

private:
    bool isExceptionPair(int particle1, int particle2) const {
        return exceptionPairs_.find(orderedPair(particle1, particle2)) != exceptionPairs_.end();
    }

    double calculate(ContextImpl& context,
                     const std::vector<Vec3>& positions,
                     std::vector<Vec3>& forces,
                     const Vec3 boxVectors[3],
                     double& dEdLambda) {
        const double lambda = getContextParameter(context, owner_->getLambdaParameterName(), 0.0);
        buildParticleStates(*owner_, lambda, particleState_);
        buildExceptionStates(*owner_, lambda, exceptionState_);

        const bool periodic = owner_->usesPeriodicBoundaryConditions();
        const bool usePME = owner_->getNonbondedMethod() == AlchemicalPMEForce::PME;
        const bool useReactionField = owner_->getNonbondedMethod() == AlchemicalPMEForce::CutoffPeriodic;
        const double cutoff = owner_->getCutoffDistance();
        const double cutoff2 = cutoff*cutoff;
        double reactionFieldK = 0.0;
        double reactionFieldC = 0.0;
        if (useReactionField) {
            const double dielectric = owner_->getReactionFieldDielectric();
            reactionFieldK = (dielectric - 1.0)/(2.0*dielectric + 1.0)/(cutoff*cutoff*cutoff);
            reactionFieldC = (3.0*dielectric)/(2.0*dielectric + 1.0)/cutoff;
        }
        double alpha = 0.0;
        int nx = 0, ny = 0, nz = 0;
        owner_->getPMEParameters(alpha, nx, ny, nz);
        if (usePME && (alpha == 0.0 || nx <= 0 || ny <= 0 || nz <= 0)) {
            throw OpenMMException("AlchemicalPMEForce PME mode requires explicit non-zero PME parameters");
        }

        double energy = 0.0;
        dEdLambda = 0.0;

        for (int i = 0; i < static_cast<int>(particleState_.size()); ++i) {
            for (int j = i + 1; j < static_cast<int>(particleState_.size()); ++j) {
                if (isExceptionPair(i, j)) {
                    continue;
                }
                Vec3 delta = positions[i] - positions[j];
                if (periodic) {
                    delta = minimumImage(delta, boxVectors);
                }
                const double r2 = dot(delta, delta);
                if ((owner_->getNonbondedMethod() != AlchemicalPMEForce::NoCutoff) && r2 >= cutoff2) {
                    continue;
                }
                const double r = std::sqrt(r2);
                if (r == 0.0) {
                    throw OpenMMException("AlchemicalPMEForce encountered two particles at identical positions");
                }
                const double invR = 1.0/r;
                const double invR2 = invR*invR;
                const auto& pi = particleState_[i];
                const auto& pj = particleState_[j];

                const double chargeProd = pi.charge*pj.charge;
                const double dChargeProd = pi.dCharge*pj.charge + pi.charge*pj.dCharge;
                if (chargeProd != 0.0 || dChargeProd != 0.0) {
                    double screening = 1.0;
                    double forceScale = invR2;
                    if (usePME) {
                        const double ar = alpha*r;
                        const double erfcTerm = std::erfc(ar);
                        screening = erfcTerm;
                        forceScale = erfcTerm*invR2 + (2.0*alpha/SqrtPi)*std::exp(-ar*ar)*invR;
                    }
                    else if (useReactionField) {
                        screening = 1.0 + reactionFieldK*r2*r - reactionFieldC*r;
                        forceScale = invR2 - 2.0*reactionFieldK*r;
                    }
                    energy += CoulombConstant*chargeProd*screening*invR;
                    dEdLambda += CoulombConstant*dChargeProd*screening*invR;
                    const Vec3 force = delta*(CoulombConstant*chargeProd*forceScale*invR);
                    forces[i] += force;
                    forces[j] -= force;
                }

                const double sigma = 0.5*(pi.sigma + pj.sigma);
                const double dSigma = 0.5*(pi.dSigma + pj.dSigma);
                const double epsilon = std::sqrt(pi.epsilon*pj.epsilon);
                double dEpsilon = 0.0;
                if (epsilon != 0.0) {
                    dEpsilon = 0.5*epsilon*(pi.dEpsilon/pi.epsilon + pj.dEpsilon/pj.epsilon);
                }
                if (sigma != 0.0 && epsilon != 0.0) {
                    const double sr = sigma*invR;
                    const double sr2 = sr*sr;
                    const double sr6 = sr2*sr2*sr2;
                    const double sr12 = sr6*sr6;
                    energy += 4.0*epsilon*(sr12 - sr6);
                    dEdLambda += 4.0*dEpsilon*(sr12 - sr6) +
                                  4.0*epsilon*(12.0*sr12 - 6.0*sr6)*(dSigma/sigma);
                    const Vec3 force = delta*(24.0*epsilon*(2.0*sr12 - sr6)*invR2);
                    forces[i] += force;
                    forces[j] -= force;
                }
            }
        }

        for (const auto& e : exceptionState_) {
            Vec3 delta = positions[e.particle1] - positions[e.particle2];
            if (periodic) {
                delta = minimumImage(delta, boxVectors);
            }
            const double r2 = dot(delta, delta);
            const double r = std::sqrt(r2);
            if (r == 0.0) {
                throw OpenMMException("AlchemicalPMEForce encountered an exception with identical particle positions");
            }
            const double invR = 1.0/r;
            const double invR2 = invR*invR;
            const auto& pi = particleState_[e.particle1];
            const auto& pj = particleState_[e.particle2];
            const double normalChargeProd = pi.charge*pj.charge;
            const double normalDChargeProd = pi.dCharge*pj.charge + pi.charge*pj.dCharge;
            if (e.chargeProd != 0.0 || e.dChargeProd != 0.0 ||
                (usePME && (normalChargeProd != 0.0 || normalDChargeProd != 0.0))) {
                double effectiveChargeProd = e.chargeProd;
                double effectiveDChargeProd = e.dChargeProd;
                double forceScale = invR2*invR;
                if (usePME) {
                    const double ar = alpha*r;
                    const double erfTerm = std::erf(ar);
                    const double derfDr = (2.0*alpha/SqrtPi)*std::exp(-ar*ar);
                    effectiveChargeProd -= normalChargeProd*erfTerm;
                    effectiveDChargeProd -= normalDChargeProd*erfTerm;
                    forceScale = effectiveChargeProd*invR2*invR + normalChargeProd*derfDr*invR2;
                }
                energy += CoulombConstant*effectiveChargeProd*invR;
                dEdLambda += CoulombConstant*effectiveDChargeProd*invR;
                const Vec3 force = usePME ?
                    delta*(CoulombConstant*forceScale) :
                    delta*(CoulombConstant*effectiveChargeProd*invR2*invR);
                forces[e.particle1] += force;
                forces[e.particle2] -= force;
            }
            if (e.sigma != 0.0 && e.epsilon != 0.0) {
                const double sr = e.sigma*invR;
                const double sr2 = sr*sr;
                const double sr6 = sr2*sr2*sr2;
                const double sr12 = sr6*sr6;
                energy += 4.0*e.epsilon*(sr12 - sr6);
                dEdLambda += 4.0*e.dEpsilon*(sr12 - sr6) +
                              4.0*e.epsilon*(12.0*sr12 - 6.0*sr6)*(e.dSigma/e.sigma);
                const Vec3 force = delta*(24.0*e.epsilon*(2.0*sr12 - sr6)*invR2);
                forces[e.particle1] += force;
                forces[e.particle2] -= force;
            }
        }

        if (usePME) {
            std::vector<double> charges(particleState_.size());
            std::vector<double> dCharges(particleState_.size());
            std::vector<int> chargeIndices;
            for (int i = 0; i < static_cast<int>(particleState_.size()); ++i) {
                charges[i] = particleState_[i].charge;
                dCharges[i] = particleState_[i].dCharge;
                if (dCharges[i] != 0.0) {
                    chargeIndices.push_back(i);
                }
            }

            const int grid[3] = {nx, ny, nz};
            double reciprocalEnergy = 0.0;
            std::vector<Vec3> reciprocalForces(particleState_.size(), Vec3());
            ReferencePME pme(alpha, static_cast<int>(particleState_.size()), grid,
                             owner_->getPMEInterpolationOrder(), 1.0);
            pme.exec(positions, reciprocalForces, charges, boxVectors, reciprocalEnergy);
            energy += reciprocalEnergy;
            for (int i = 0; i < static_cast<int>(reciprocalForces.size()); ++i) {
                forces[i] += reciprocalForces[i];
            }

            if (!chargeIndices.empty()) {
                std::vector<double> chargeDerivatives(chargeIndices.size(), 0.0);
                pme.exec_charge_derivatives(positions, chargeDerivatives, chargeIndices, charges, boxVectors);
                for (int i = 0; i < static_cast<int>(chargeIndices.size()); ++i) {
                    dEdLambda += chargeDerivatives[i]*dCharges[chargeIndices[i]];
                }
            }

            double selfEnergy = 0.0;
            double selfDerivative = 0.0;
            double totalCharge = 0.0;
            double totalDCharge = 0.0;
            for (int i = 0; i < static_cast<int>(particleState_.size()); ++i) {
                selfEnergy += charges[i]*charges[i];
                selfDerivative += 2.0*charges[i]*dCharges[i];
                totalCharge += charges[i];
                totalDCharge += dCharges[i];
            }
            energy -= CoulombConstant*alpha*selfEnergy/SqrtPi;
            dEdLambda -= CoulombConstant*alpha*selfDerivative/SqrtPi;

            const double volume = dot(boxVectors[0], cross(boxVectors[1], boxVectors[2]));
            if (volume <= 0.0) {
                throw OpenMMException("AlchemicalPMEForce requires a positive periodic box volume");
            }
            energy -= CoulombConstant*Pi*totalCharge*totalCharge/(2.0*volume*alpha*alpha);
            dEdLambda -= CoulombConstant*Pi*totalCharge*totalDCharge/(volume*alpha*alpha);
        }

        return energy;
    }

    const AlchemicalPMEForce* owner_;
    std::vector<Vec3> localForces_;
    std::vector<ParticleState> particleState_;
    std::vector<ExceptionState> exceptionState_;
    std::set<std::pair<int, int>> exceptionPairs_;
};

class ReferenceAlchemicalPMEKernelFactory final : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const override {
        if (name == CalcAlchemicalPMEForceKernel::Name()) {
            return new ReferenceCalcAlchemicalPMEForceKernel(name, platform);
        }
        throw OpenMMException("Tried to create an illegal AlchemicalPMEForce kernel");
    }
};

std::mutex registerFactoriesMutex;

} // namespace

void OpenMM::registerAlchemicalPMEKernelFactories()
{
    std::lock_guard<std::mutex> lock(registerFactoriesMutex);
    for (int i = 0; i < Platform::getNumPlatforms(); ++i) {
        Platform& platform = Platform::getPlatform(i);
        if ((platform.getName() == "Reference" || platform.getName() == "CPU") &&
            !platform.supportsKernels({CalcAlchemicalPMEForceKernel::Name()})) {
            platform.registerKernelFactory(CalcAlchemicalPMEForceKernel::Name(),
                                           new ReferenceAlchemicalPMEKernelFactory());
        }
        if (platform.getName() == "CUDA" &&
            !platform.supportsKernels({CalcAlchemicalPMEForceKernel::Name()})) {
            registerCudaAlchemicalPMEKernelFactory(platform);
        }
    }
}
