#ifndef OPENMM_INTERNAL_COLVARS_FORCE_IMPL_H
#define OPENMM_INTERNAL_COLVARS_FORCE_IMPL_H

#include "openmm/ColvarsForce.h"
#include "openmm/internal/CustomCPPForceImpl.h"

#include <memory>
#include <map>
#include <vector>

namespace OpenMM {

class OpenMMColvarProxy;

class ColvarsForceImpl final : public CustomCPPForceImpl {
public:
    explicit ColvarsForceImpl(const ColvarsForce& owner);

    const ColvarsForce& getOwner() const;

    void initialize(ContextImpl& context) override;
    double calcForcesAndEnergy(ContextImpl& context,
                               bool includeForces,
                               bool includeEnergy,
                               int groups) override;
    double computeForce(ContextImpl& context,
                        const std::vector<Vec3>& positions,
                        std::vector<Vec3>& forces) override;
    std::map<std::string, double> getDefaultParameters() override;

private:
    const ColvarsForce& owner_;
    std::unique_ptr<OpenMMColvarProxy> proxy_;
    bool usesDeferredCudaComputation_;
};

} // namespace OpenMM

#endif
