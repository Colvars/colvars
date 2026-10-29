#include "ColvarsForceImpl.h"
#include "colvarproxy_openmm.h"

#include "openmm/CustomExternalForce.h"
#include "openmm/internal/ContextImpl.h"

#include <algorithm>
#include <cstdint>

#ifdef OPENMM_COLVARS_ENABLE_CUDA_KERNELS
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/cuda/CudaPlatform.h"
#endif

using namespace OpenMM;
using namespace ColvarsPlugin;

#ifdef OPENMM_COLVARS_ENABLE_CUDA_KERNELS
namespace {

const char* deferredColvarsForceSource = R"(
KERNEL void addForces(GLOBAL const real* RESTRICT forces,
                      GLOBAL mm_long* RESTRICT forceBuffers,
                      GLOBAL int* RESTRICT atomIndex) {
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        int index = atomIndex[atom];
        forceBuffers[atom] += (mm_long) (forces[3*index]*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] +=
            (mm_long) (forces[3*index+1]*0x100000000);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] +=
            (mm_long) (forces[3*index+2]*0x100000000);
    }
}
)";

class DeferredCudaColvarsComputation final : public ComputeContext::ForcePostComputation {
public:
    DeferredCudaColvarsComputation(ContextImpl& context,
                                   CudaContext& computeContext,
                                   ColvarsForceImpl& force)
        : context_(context),
          computeContext_(computeContext),
          cudaContext_(computeContext),
          force_(force),
          forceGroupFlag_(std::uint32_t{1} <<
                          static_cast<unsigned>(force.getOwner().getForceGroup()))
    {
        ContextSelector selector(computeContext_);
        const int numParticles = context_.getSystem().getNumParticles();
        const int paddedNumParticles = computeContext_.getPaddedNumAtoms();
        positions_.resize(numParticles);
        forces_.resize(numParticles);
        floatForces_.resize(3 * numParticles);
        const int elementSize =
            computeContext_.getUseDoublePrecision() ? sizeof(double) : sizeof(float);
        forceArray_.initialize(computeContext_,
                               3 * numParticles,
                               elementSize,
                               "deferredColvarsForces");
        std::map<std::string, std::string> defines;
        defines["NUM_ATOMS"] = computeContext_.intToString(numParticles);
        defines["PADDED_NUM_ATOMS"] =
            computeContext_.intToString(computeContext_.getPaddedNumAtoms());
        ComputeProgram program =
            computeContext_.compileProgram(deferredColvarsForceSource, defines);
        addForcesKernel_ = program->createKernel("addForces");
        addForcesKernel_->addArg(forceArray_);
        addForcesKernel_->addArg(computeContext_.getLongForceBuffer());
        addForcesKernel_->addArg(computeContext_.getAtomIndexArray());
    }

    double computeForceAndEnergy(bool includeForces,
                                 bool includeEnergy,
                                 int groups) override
    {
        (void) includeEnergy;
        if ((static_cast<std::uint32_t>(groups) & forceGroupFlag_) == 0) {
            return 0.0;
        }

        cudaContext_.flushQueue();
        context_.getPositions(positions_);
        const double energy = force_.computeForce(context_, positions_, forces_);
        if (!includeForces) {
            return energy;
        }

        ContextSelector selector(computeContext_);
        if (computeContext_.getUseDoublePrecision()) {
            forceArray_.upload(reinterpret_cast<const double*>(forces_.data()));
        } else {
            for (int i = 0; i < static_cast<int>(forces_.size()); ++i) {
                floatForces_[3 * i] = static_cast<float>(forces_[i][0]);
                floatForces_[3 * i + 1] = static_cast<float>(forces_[i][1]);
                floatForces_[3 * i + 2] = static_cast<float>(forces_[i][2]);
            }
            forceArray_.upload(floatForces_);
        }
        addForcesKernel_->execute(computeContext_.getNumAtoms());
        return energy;
    }

private:
    ContextImpl& context_;
    ComputeContext& computeContext_;
    CudaContext& cudaContext_;
    ColvarsForceImpl& force_;
    std::uint32_t forceGroupFlag_;
    ComputeArray forceArray_;
    ComputeKernel addForcesKernel_;
    std::vector<Vec3> positions_;
    std::vector<Vec3> forces_;
    std::vector<float> floatForces_;
};

class KeepDeferredColvarsLast final : public ComputeContext::ForcePreComputation {
public:
    KeepDeferredColvarsLast(ComputeContext& context,
                            ComputeContext::ForcePostComputation* colvars)
        : context_(context), colvars_(colvars)
    {
    }

    void computeForceAndEnergy(bool includeForces,
                               bool includeEnergy,
                               int groups) override
    {
        (void) includeForces;
        (void) includeEnergy;
        (void) groups;
        auto& computations = context_.getPostComputations();
        const auto found = std::find(computations.begin(), computations.end(), colvars_);
        if (found != computations.end() && found + 1 != computations.end()) {
            computations.erase(found);
            computations.push_back(colvars_);
        }
    }

private:
    ComputeContext& context_;
    ComputeContext::ForcePostComputation* colvars_;
};

} // namespace
#endif

ColvarsForceImpl::ColvarsForceImpl(const ColvarsForce& owner)
    : CustomCPPForceImpl(owner),
      owner_(owner),
      usesDeferredCudaComputation_(false)
{
}

namespace {

bool isOwnerLastForce(const System& system, const Force& owner)
{
    const int numForces = system.getNumForces();
    return numForces > 0 && &system.getForce(numForces - 1) == &owner;
}

bool isSafeAfterColvarsForTotalForces(const Force& force)
{
    const auto* external = dynamic_cast<const CustomExternalForce*>(&force);
    return external != nullptr && external->getEnergyFunction() == "0";
}

bool hasOnlySafeForcesAfterOwner(const System& system, const Force& owner)
{
    bool afterOwner = false;
    for (int forceIndex = 0; forceIndex < system.getNumForces(); ++forceIndex) {
        const Force& force = system.getForce(forceIndex);
        if (&force == &owner) {
            afterOwner = true;
            continue;
        }
        if (afterOwner && !isSafeAfterColvarsForTotalForces(force)) {
            return false;
        }
    }
    return afterOwner;
}

void validateTotalForceOrdering(const ContextImpl& context,
                                const ColvarsForce& owner,
                                const colvarproxy_openmm& proxy)
{
    if (!proxy.total_forces_enabled() || isOwnerLastForce(context.getSystem(), owner) ||
        hasOnlySafeForcesAfterOwner(context.getSystem(), owner)) {
        return;
    }
    throw OpenMMException(
        "Colvars configurations that request total forces, such as ABF/eABF, "
        "require ColvarsForce to execute after every physical Force on Reference "
        "and CPU platforms so the current OpenMM force buffer is complete.");
}

void validateTotalForceConfiguration(ContextImpl& context,
                                     const ColvarsForce& owner,
                                     const colvarproxy_openmm& proxy,
                                     bool usesDeferredCudaComputation)
{
    if (!proxy.total_forces_enabled()) {
        return;
    }

    int colvarsForces = 0;
    for (int forceIndex = 0; forceIndex < context.getSystem().getNumForces(); ++forceIndex) {
        if (dynamic_cast<const ColvarsForce*>(&context.getSystem().getForce(forceIndex)) != nullptr) {
            ++colvarsForces;
        }
    }
    if (colvarsForces > 1) {
        throw OpenMMException(
            "A System containing a Colvars configuration that requests total forces "
            "may contain only one ColvarsForce. Put multiple Colvars biases in one configuration.");
    }

    const std::string platformName = context.getPlatform().getName();
    if (usesDeferredCudaComputation) {
        return;
    }
    if (platformName != "Reference" && platformName != "CPU") {
        throw OpenMMException(
            "Colvars total-force sampling is unsupported on platform \"" + platformName +
            "\" because the completed force buffer is unavailable. Use Reference or CPU, "
            "or build with OPENMM_COLVARS_ENABLE_CUDA_KERNELS=ON and use CUDA.");
    }
    validateTotalForceOrdering(context, owner, proxy);
}

void validateNonDeferredForceGroup(const ColvarsForce& owner)
{
    if (owner.getForceGroup() == 31) {
        throw OpenMMException(
            "ColvarsForce group 31 is unsupported on the non-deferred CustomCPPForce path "
            "because this OpenMM dependency uses a signed force-group shift. Use group 0-30, "
            "or use a deferred CUDA build.");
    }
}

} // namespace

const ColvarsForce& ColvarsForceImpl::getOwner() const
{
    return owner_;
}

void ColvarsForceImpl::initialize(ContextImpl& context)
{
#ifdef OPENMM_COLVARS_ENABLE_CUDA_KERNELS
    if (context.getPlatform().getName() == "CUDA") {
        auto* data =
            reinterpret_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
        if (data == nullptr || data->contexts.size() != 1) {
            throw OpenMMException(
                "CUDA ColvarsForce currently supports a single CUDA context only");
        }
        proxy_ = std::make_unique<colvarproxy_openmm>(owner_, context.getSystem());
        validateTotalForceConfiguration(context, owner_, *proxy_, true);
        CudaContext& computeContext = *data->contexts[0];
        auto* deferred =
            new DeferredCudaColvarsComputation(context, computeContext, *this);
        computeContext.addPostComputation(deferred);
        computeContext.addPreComputation(
            new KeepDeferredColvarsLast(computeContext, deferred));
        usesDeferredCudaComputation_ = true;
        return;
    }
#endif
    validateNonDeferredForceGroup(owner_);
    CustomCPPForceImpl::initialize(context);
    proxy_ = std::make_unique<colvarproxy_openmm>(owner_, context.getSystem());
    validateTotalForceConfiguration(context, owner_, *proxy_, false);
}

double ColvarsForceImpl::calcForcesAndEnergy(ContextImpl& context,
                                             bool includeForces,
                                             bool includeEnergy,
                                             int groups)
{
    if (usesDeferredCudaComputation_) {
        return 0.0;
    }
    return CustomCPPForceImpl::calcForcesAndEnergy(
        context, includeForces, includeEnergy, groups);
}

double ColvarsForceImpl::computeForce(ContextImpl& context,
                                      const std::vector<Vec3>& positions,
                                      std::vector<Vec3>& forces)
{
    std::fill(forces.begin(), forces.end(), Vec3());
    if (proxy_ == nullptr) {
        proxy_ = std::make_unique<colvarproxy_openmm>(owner_, context.getSystem());
    }

    const double energy = proxy_->calculate(context, positions);
    proxy_->addForcesToOpenMM(forces);
    return energy;
}

std::map<std::string, double> ColvarsForceImpl::getDefaultParameters()
{
    std::map<std::string, double> parameters;
    if (owner_.getAlchemyEnabled()) {
        parameters[owner_.getLambdaParameterName()] = 0.0;
        parameters[owner_.getLambdaDerivativeParameterName()] = 0.0;
    }
    return parameters;
}
