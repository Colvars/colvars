// Adapted from Colvars PR 979. See ../UPSTREAM.md and ../LICENSE.
#include "ColvarsForceImpl.h"
#include "colvarproxy_openmm.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include <algorithm>
#include <cstdint>
using namespace OpenMM;
using namespace ColvarsPlugin;

ColvarsForceImpl::ColvarsForceImpl(const ColvarsForce& owner)
    : CustomCPPForceImpl(owner), owner_(owner) {}
ColvarsForceImpl::~ColvarsForceImpl() {
    OpenMMColvarsLock lock;
    cuda_.reset();
    proxy_.reset();
}
void ColvarsForceImpl::initialize(ContextImpl& context) {
    OpenMMColvarsLock lock;
    platform_ = context.getPlatform().getName();
    if (platform_ == "CUDA") {
#ifdef OPENMM_COLVARS_ENABLE_CUDA
        cuda_ = createCudaColvarsKernel(context, owner_);
        return;
#else
        throw OpenMMException("Colvars was built without CUDA support; rebuild with OPENMM_COLVARS_ENABLE_CUDA=ON");
#endif
    }
    if (platform_ != "Reference" && platform_ != "CPU")
        throw OpenMMException("Colvars supports OpenMM Reference, CPU, and CUDA platforms");
    // Older OpenMM CustomCPP kernels use a signed shift for the group mask.
    if (owner_.getForceGroup() == 31)
        throw OpenMMException("Colvars force group 31 requires CUDA; use a group from 0 through 30 on CPU/Reference");
    CustomCPPForceImpl::initialize(context);
    proxy_ = std::make_unique<colvarproxy_openmm>(owner_, context.getSystem());
    proxy_->initializeColvars();
}
double ColvarsForceImpl::calcForcesAndEnergy(ContextImpl& c, bool forces, bool energy, int groups) {
    OpenMMColvarsLock lock(&proxy());
    if (!(static_cast<std::uint32_t>(groups) & (std::uint32_t(1) << owner_.getForceGroup()))) return 0.0;
    if (!forces && !energy) return 0.0;
    if (cuda_) return cuda_->execute(c, forces, energy);
    return CustomCPPForceImpl::calcForcesAndEnergy(c, forces, energy, groups);
}
double ColvarsForceImpl::computeForce(ContextImpl& c, const std::vector<Vec3>& positions,
                                     std::vector<Vec3>& forces) {
    OpenMMColvarsLock lock(proxy_.get());
    std::fill(forces.begin(), forces.end(), Vec3());
    const double energy = proxy_->calculate(c, positions);
    proxy_->addForcesToOpenMM(forces);
    return energy;
}
colvarproxy_openmm& ColvarsForceImpl::proxy() { return cuda_ ? cuda_->proxy() : *proxy_; }
std::string ColvarsForceImpl::backend() const {
    return cuda_ ? cuda_->backend() : platform_ + "/Colvars CPU";
}
