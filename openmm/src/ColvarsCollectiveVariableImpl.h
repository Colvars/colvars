#ifndef OPENMM_COLVARS_COLLECTIVE_VARIABLE_IMPL_H
#define OPENMM_COLVARS_COLLECTIVE_VARIABLE_IMPL_H

#include "ColvarsCollectiveVariable.h"
#include "ColvarsForce.h"
#include "openmm/internal/CustomCPPForceImpl.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

class colvarproxy_openmm;

namespace ColvarsPlugin {

class ColvarsCollectiveVariableImpl final : public OpenMM::CustomCPPForceImpl {
public:
    explicit ColvarsCollectiveVariableImpl(const ColvarsCollectiveVariable& owner);
    ~ColvarsCollectiveVariableImpl() override;

    const OpenMM::Force& getOwner() const override;
    void initialize(OpenMM::ContextImpl& context) override;
    double computeForce(OpenMM::ContextImpl& context,
                        const std::vector<OpenMM::Vec3>& positions,
                        std::vector<OpenMM::Vec3>& forces) override;
    std::map<std::string, double> getDefaultParameters() override;

private:
    const ColvarsCollectiveVariable& owner_;
    // Holds the Colvars settings that colvarproxy_openmm reads; it is never
    // added to a System
    ColvarsForce settings_;
    std::unique_ptr<colvarproxy_openmm> proxy_;
};

} // namespace ColvarsPlugin

#endif
