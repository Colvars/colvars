#ifndef OPENMM_COLVARS_FORCE_IMPL_H
#define OPENMM_COLVARS_FORCE_IMPL_H
#include "ColvarsForce.h"
#include "CudaColvarsKernel.h"
#include "openmm/internal/CustomCPPForceImpl.h"
#include <memory>
class colvarproxy_openmm;
namespace ColvarsPlugin {
class ColvarsForceImpl : public OpenMM::CustomCPPForceImpl {
public:
    explicit ColvarsForceImpl(const ColvarsForce&);
    ~ColvarsForceImpl() override;
    const ColvarsForce& getOwner() const override { return owner_; }
    void initialize(OpenMM::ContextImpl&) override;
    double calcForcesAndEnergy(OpenMM::ContextImpl&, bool, bool, int) override;
    double computeForce(OpenMM::ContextImpl&, const std::vector<OpenMM::Vec3>&,
                         std::vector<OpenMM::Vec3>&) override;
    colvarproxy_openmm& proxy();
    std::string backend() const;
private:
    const ColvarsForce& owner_;
    std::unique_ptr<colvarproxy_openmm> proxy_;
    std::unique_ptr<CudaColvarsKernel> cuda_;
    std::string platform_;
};
}
#endif
