#ifndef OPENMM_ALCHEMICAL_PME_FORCE_IMPL_H
#define OPENMM_ALCHEMICAL_PME_FORCE_IMPL_H

#include "AlchemicalPMEForce.h"
#include "openmm/Kernel.h"
#include "openmm/internal/ForceImpl.h"

#include <map>
#include <vector>

namespace ColvarsPlugin {

class AlchemicalPMEForceImpl final : public OpenMM::ForceImpl {
public:
    explicit AlchemicalPMEForceImpl(const AlchemicalPMEForce& owner);

    const AlchemicalPMEForce& getOwner() const;
    void initialize(OpenMM::ContextImpl& context) override;
    void updateContextState(OpenMM::ContextImpl& context, bool& forcesInvalid) override;
    double calcForcesAndEnergy(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy, int groups) override;
    std::map<std::string, double> getDefaultParameters() override;
    std::vector<std::string> getKernelNames() override;

private:
    const AlchemicalPMEForce& owner_;
    OpenMM::Kernel kernel_;
};

} // namespace ColvarsPlugin

#endif
