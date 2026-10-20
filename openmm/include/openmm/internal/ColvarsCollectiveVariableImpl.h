#ifndef OPENMM_INTERNAL_COLVARS_COLLECTIVE_VARIABLE_IMPL_H
#define OPENMM_INTERNAL_COLVARS_COLLECTIVE_VARIABLE_IMPL_H

#include "openmm/ColvarsCollectiveVariable.h"
#include "openmm/ColvarsForce.h"
#include "openmm/internal/CustomCPPForceImpl.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace OpenMM {

class OpenMMColvarProxy;

class ColvarsCollectiveVariableImpl final : public CustomCPPForceImpl {
public:
    explicit ColvarsCollectiveVariableImpl(const ColvarsCollectiveVariable& owner);
    ~ColvarsCollectiveVariableImpl() override;

    const Force& getOwner() const override;
    void initialize(ContextImpl& context) override;
    double computeForce(ContextImpl& context,
                        const std::vector<Vec3>& positions,
                        std::vector<Vec3>& forces) override;
    std::map<std::string, double> getDefaultParameters() override;

private:
    const ColvarsCollectiveVariable& owner_;
    // Holds the Colvars settings that OpenMMColvarProxy reads; it is never
    // added to a System
    ColvarsForce settings_;
    std::unique_ptr<OpenMMColvarProxy> proxy_;
};

} // namespace OpenMM

#endif
