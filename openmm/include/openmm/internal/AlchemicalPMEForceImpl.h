#ifndef OPENMM_INTERNAL_ALCHEMICAL_PME_FORCE_IMPL_H
#define OPENMM_INTERNAL_ALCHEMICAL_PME_FORCE_IMPL_H

#include "openmm/AlchemicalPMEForce.h"
#include "openmm/Kernel.h"
#include "openmm/internal/ForceImpl.h"

#include <map>
#include <vector>

namespace OpenMM {

class AlchemicalPMEForceImpl final : public ForceImpl {
public:
    explicit AlchemicalPMEForceImpl(const AlchemicalPMEForce& owner);

    const AlchemicalPMEForce& getOwner() const;
    void initialize(ContextImpl& context) override;
    void updateContextState(ContextImpl& context, bool& forcesInvalid) override;
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) override;
    std::map<std::string, double> getDefaultParameters() override;
    std::vector<std::string> getKernelNames() override;

private:
    const AlchemicalPMEForce& owner_;
    Kernel kernel_;
};

} // namespace OpenMM

#endif
