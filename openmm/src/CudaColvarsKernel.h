#ifndef OPENMM_CUDA_COLVARS_KERNEL_H
#define OPENMM_CUDA_COLVARS_KERNEL_H
#include <memory>
#include <string>
class colvarproxy_openmm;
namespace OpenMM { class ContextImpl; }
namespace ColvarsPlugin {
class ColvarsForce;
class CudaColvarsKernel {
public:
    virtual ~CudaColvarsKernel() = default;
    virtual double execute(OpenMM::ContextImpl&, bool includeForces, bool includeEnergy) = 0;
    virtual colvarproxy_openmm& proxy() = 0;
    virtual std::string backend() const = 0;
};
std::unique_ptr<CudaColvarsKernel> createCudaColvarsKernel(
    OpenMM::ContextImpl&, const ColvarsForce&);
}
#endif
