#ifndef OPENMM_COLVARS_FORCE_IMPL_H
#define OPENMM_COLVARS_FORCE_IMPL_H

#include "ColvarsForce.h"
#include "openmm/internal/CustomCPPForceImpl.h"

#include <memory>
#include <map>
#include <vector>

class colvarproxy_openmm;

namespace ColvarsPlugin {

class ColvarsForceImpl final : public OpenMM::CustomCPPForceImpl {
public:
    explicit ColvarsForceImpl(const ColvarsForce& owner);

    const ColvarsForce& getOwner() const;

    void initialize(OpenMM::ContextImpl& context) override;
    double calcForcesAndEnergy(OpenMM::ContextImpl& context,
                               bool includeForces,
                               bool includeEnergy,
                               int groups) override;
    double computeForce(OpenMM::ContextImpl& context,
                        const std::vector<OpenMM::Vec3>& positions,
                        std::vector<OpenMM::Vec3>& forces) override;
    std::map<std::string, double> getDefaultParameters() override;

private:
    const ColvarsForce& owner_;
    std::unique_ptr<colvarproxy_openmm> proxy_;
    bool usesDeferredCudaComputation_;
};

} // namespace ColvarsPlugin

#endif
