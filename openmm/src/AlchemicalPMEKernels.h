#ifndef OPENMM_ALCHEMICAL_PME_KERNELS_H
#define OPENMM_ALCHEMICAL_PME_KERNELS_H

#include "AlchemicalPMEForce.h"
#include "openmm/KernelImpl.h"

namespace OpenMM {
class ContextImpl;
}

namespace ColvarsPlugin {

class CalcAlchemicalPMEForceKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "ColvarsCalcAlchemicalPMEForce";
    }

    CalcAlchemicalPMEForceKernel(std::string name, const OpenMM::Platform& platform)
        : OpenMM::KernelImpl(name, platform) {
    }

    virtual void initialize(OpenMM::ContextImpl& context, const AlchemicalPMEForce& force) = 0;
    virtual double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
};

void registerAlchemicalPMEKernelFactories();
#ifdef OPENMM_COLVARS_ENABLE_CUDA_KERNELS
void registerCudaAlchemicalPMEKernelFactory(OpenMM::Platform& platform);
#endif

} // namespace ColvarsPlugin

#endif
