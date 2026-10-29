#ifndef OPENMM_ALCHEMICAL_PME_KERNELS_H
#define OPENMM_ALCHEMICAL_PME_KERNELS_H

#include "AlchemicalPMEForce.h"
#include "openmm/KernelImpl.h"

namespace OpenMM {

class ContextImpl;

class CalcAlchemicalPMEForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcAlchemicalPMEForce";
    }

    CalcAlchemicalPMEForceKernel(std::string name, const Platform& platform)
        : KernelImpl(name, platform) {
    }

    virtual void initialize(ContextImpl& context, const AlchemicalPMEForce& force) = 0;
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
};

void registerAlchemicalPMEKernelFactories();
#ifdef OPENMM_COLVARS_ENABLE_CUDA_KERNELS
void registerCudaAlchemicalPMEKernelFactory(Platform& platform);
#endif

} // namespace OpenMM

#endif
