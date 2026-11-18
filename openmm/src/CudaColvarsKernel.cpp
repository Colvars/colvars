// OpenMM CUDA adapter for the Colvars GPU implementation in ../src/cuda.
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "CudaColvarsKernel.h"
#include "ColvarsForce.h"
#include "colvarproxy_openmm.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/ComputeKernel.h"
#include "openmm/common/ComputeProgram.h"
#include "openmm/cuda/CudaArray.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/cuda/CudaPlatform.h"
#include <algorithm>
#include <cstdint>
#include <memory>
#include <sstream>
#include <vector>

using namespace OpenMM;
using namespace ColvarsPlugin;

namespace {

void checkCuda(cudaError_t result, const char* operation) {
    if (result != cudaSuccess)
        throw OpenMMException(std::string("Colvars CUDA ")+operation+": "+cudaGetErrorString(result));
}

// OpenMM keeps positions in reordered real4 arrays and, in mixed precision,
// a separate float4 correction.  Colvars expects selected atoms in double SoA.
// Cell offsets undo OpenMM's internal molecule wrapping, exactly as getPositions
// does.  The Colvars components then apply their own minimum-image convention.
const char* adapterSource = R"CUDA(
KERNEL void gatherColvarsPositions(GLOBAL const real4* RESTRICT posq,
                                  GLOBAL const real4* RESTRICT correction,
                                  GLOBAL const int* RESTRICT selected,
                                  GLOBAL const int4* RESTRICT cellOffsets,
                                  GLOBAL double* RESTRICT positions,
                                  double4 boxA, double4 boxB, double4 boxC) {
    for (int i = GLOBAL_ID; i < NUM_SELECTED; i += GLOBAL_SIZE) {
        const int atom = selected[i];
        const real4 p = posq[atom];
        double x = (double) p.x, y = (double) p.y, z = (double) p.z;
#ifdef USE_MIXED_PRECISION
        const real4 c = correction[atom];
        x += (double) c.x;
        y += (double) c.y;
        z += (double) c.z;
#endif
        const int4 offset = cellOffsets[i];
        positions[i] = x-boxA.x*offset.x-boxB.x*offset.y-boxC.x*offset.z;
        positions[i+NUM_SELECTED] = y-boxA.y*offset.x-boxB.y*offset.y-boxC.y*offset.z;
        positions[i+2*NUM_SELECTED] = z-boxA.z*offset.x-boxB.z*offset.y-boxC.z*offset.z;
    }
}

KERNEL void addColvarsForces(GLOBAL const double* RESTRICT forces,
                            GLOBAL const int* RESTRICT selected,
                            GLOBAL mm_ulong* RESTRICT forceBuffer) {
    for (int i = GLOBAL_ID; i < NUM_SELECTED; i += GLOBAL_SIZE) {
        const int atom = selected[i];
        // Both packages use kJ/mol/nm.  OpenMM's accumulator has 32 fractional
        // bits.  Atomics also make this safe if an OpenMM force uses another
        // stream concurrently (for example PME).
        ATOMIC_ADD(&forceBuffer[atom],
                   (mm_ulong) ((mm_long) (forces[i]*4294967296.0)));
        ATOMIC_ADD(&forceBuffer[atom+PADDED_NUM_ATOMS],
                   (mm_ulong) ((mm_long) (forces[i+NUM_SELECTED]*4294967296.0)));
        ATOMIC_ADD(&forceBuffer[atom+2*PADDED_NUM_ATOMS],
                   (mm_ulong) ((mm_long) (forces[i+2*NUM_SELECTED]*4294967296.0)));
    }
}
)CUDA";

class CudaProxy final : public colvarproxy_openmm {
public:
    CudaProxy(const ColvarsForce& owner, const System& system, CudaContext& cuda)
        : colvarproxy_openmm(owner, system), cuda_(cuda), allowGpu_(owner.getUseGpu()) {
        initializeColvars();
    }
    smp_mode_t get_preferred_smp_mode() const override {
        return allowGpu_ ? smp_mode_t::gpu : smp_mode_t::none;
    }
    std::vector<smp_mode_t> get_available_smp_modes() const override {
        if (allowGpu_) return {smp_mode_t::none, smp_mode_t::gpu};
        return {smp_mode_t::none};
    }
    int set_smp_mode(smp_mode_t mode) override {
        if (mode != smp_mode_t::none && !(allowGpu_ && mode == smp_mode_t::gpu))
            return COLVARS_NOT_IMPLEMENTED;
        smp_mode = mode;
        support_gpu = (mode == smp_mode_t::gpu);
        return COLVARS_OK;
    }
    cudaStream_t get_default_stream() override {
        return reinterpret_cast<cudaStream_t>(cuda_.getCurrentStream());
    }
    cvm::real* proxy_atoms_positions_gpu() override { return positions_; }
    cvm::real* proxy_atoms_new_colvar_forces_gpu() override { return forces_; }
    void setDeviceBuffers(double* positions, double* forces) {
        positions_ = positions;
        forces_ = forces;
        if (has_gpu_support() && module().proxy_gpu_buffers_reallocated_done() != COLVARS_OK)
            throw OpenMMException("Unable to bind Colvars GPU atom buffers");
    }
private:
    CudaContext& cuda_;
    bool allowGpu_;
    double* positions_ = nullptr;
    double* forces_ = nullptr;
};

// ComputeContext owns listeners.  Use shared state instead of a pointer to the
// kernel so destruction order cannot leave a dangling callback.
struct MappingState { bool dirty = true; };
class MappingListener final : public ComputeContext::ReorderListener {
public:
    explicit MappingListener(std::shared_ptr<MappingState> state) : state_(std::move(state)) {}
    void execute() override { state_->dirty = true; }
private:
    std::shared_ptr<MappingState> state_;
};

class PinnedBuffer {
public:
    ~PinnedBuffer() { if (data_) cudaFreeHost(data_); }
    void allocate(size_t count) {
        checkCuda(cudaHostAlloc(reinterpret_cast<void**>(&data_), count*sizeof(double), cudaHostAllocDefault),
                  "allocating pinned fallback buffer");
    }
    double* data() { return data_; }
private:
    double* data_ = nullptr;
};

class ResidentCudaColvarsKernel final : public CudaColvarsKernel {
public:
    ResidentCudaColvarsKernel(ContextImpl& context, const ColvarsForce& owner, CudaContext& cuda)
        : cuda_(cuda), mapping_(std::make_shared<MappingState>()) {
        OpenMMColvarsLock lock;
        ContextSelector selector(cuda_);
        try {
            proxy_ = std::make_unique<CudaProxy>(owner, context.getSystem(), cuda_);
            selectedCount_ = static_cast<int>(proxy_->requestedAtoms().size());
            const size_t size = std::max(1, selectedCount_);
            positions_.initialize<double>(cuda_, 3*size, "colvarsSelectedPositions");
            forces_.initialize<double>(cuda_, 3*size, "colvarsSelectedForces");
            selected_.initialize<int>(cuda_, size, "colvarsSelectedIndices");
            offsets_.initialize<mm_int4>(cuda_, size, "colvarsSelectedCellOffsets");
            reordered_.resize(size);
            cellOffsets_.resize(size, mm_int4(0, 0, 0, 0));
            inverseOrder_.resize(cuda_.getNumAtoms());
            if (!proxy_->has_gpu_support()) {
                hostPositions_.allocate(3*size);
                hostForces_.allocate(3*size);
            }
            proxy_->setDeviceBuffers(reinterpret_cast<double*>(positions_.getDevicePointer()),
                                     reinterpret_cast<double*>(forces_.getDevicePointer()));
            std::map<std::string, std::string> defines;
            defines["NUM_SELECTED"] = std::to_string(selectedCount_);
            defines["PADDED_NUM_ATOMS"] = std::to_string(cuda_.getPaddedNumAtoms());
            ComputeProgram program = cuda_.compileProgram(adapterSource, defines);
            gather_ = program->createKernel("gatherColvarsPositions");
            gather_->addArg(cuda_.getPosq());
            // Use a valid dummy pointer when the mixed-precision branch is absent.
            gather_->addArg(cuda_.getUseMixedPrecision() ? cuda_.getPosqCorrection() : cuda_.getPosq());
            gather_->addArg(selected_);
            gather_->addArg(offsets_);
            gather_->addArg(positions_);
            gather_->addArg(mm_double4(0, 0, 0, 0));
            gather_->addArg(mm_double4(0, 0, 0, 0));
            gather_->addArg(mm_double4(0, 0, 0, 0));
            scatter_ = program->createKernel("addColvarsForces");
            scatter_->addArg(forces_);
            scatter_->addArg(selected_);
            scatter_->addArg(cuda_.getLongForceBuffer());
            cuda_.addReorderListener(new MappingListener(mapping_));
            updateMapping();
            describeBackend();
            proxy_->log(backend_+"\n");
        } catch (...) {
            // Member destruction follows local-variable destruction when a
            // constructor throws.  Release the Colvars module while this
            // selector still makes its CUDA context current, since its graphs
            // and streams are not wrapped in OpenMM context selectors.
            proxy_.reset();
            throw;
        }
    }

    ~ResidentCudaColvarsKernel() override {
        // The module destroys its graphs, group streams, and pinned buffers.
        // Make the correct OpenMM CUDA context current while doing that.
        OpenMMColvarsLock lock(proxy_.get());
        ContextSelector selector(cuda_);
        cudaStreamSynchronize(reinterpret_cast<cudaStream_t>(cuda_.getCurrentStream()));
        proxy_.reset();
    }

    double execute(ContextImpl& context, bool includeForces, bool includeEnergy) override {
        if (!includeForces && !includeEnergy) return 0.0;
        OpenMMColvarsLock lock(proxy_.get());
        ContextSelector selector(cuda_);
        proxy_->beginEvaluation(context);
        updateMapping();
        const cudaStream_t stream = proxy_->get_default_stream();
        if (selectedCount_ > 0) {
            Vec3 a, b, c;
            cuda_.getPeriodicBoxVectors(a, b, c);
            gather_->setArg(5, mm_double4(a[0], a[1], a[2], 0));
            gather_->setArg(6, mm_double4(b[0], b[1], b[2], 0));
            gather_->setArg(7, mm_double4(c[0], c[1], c[2], 0));
            gather_->execute(selectedCount_);
        }
        if (proxy_->has_gpu_support()) {
            checkCuda(cudaMemsetAsync(reinterpret_cast<void*>(forces_.getDevicePointer()), 0,
                                      3*std::max(1, selectedCount_)*sizeof(double), stream), "clearing selected forces");
            checkCuda(cudaEventRecord(proxy_->get_event(colvarproxy_gpu::event_type::copy_atoms), stream),
                      "recording selected positions event");
            checkCuda(cudaEventRecord(proxy_->get_event(colvarproxy_gpu::event_type::update_lattice), stream),
                      "recording periodic box event");
        } else if (selectedCount_ > 0) {
            checkCuda(cudaMemcpyAsync(hostPositions_.data(), reinterpret_cast<void*>(positions_.getDevicePointer()),
                                      3*selectedCount_*sizeof(double), cudaMemcpyDeviceToHost, stream),
                      "downloading selected positions");
            checkCuda(cudaStreamSynchronize(stream), "waiting for selected positions");
            auto& host = *proxy_->modify_atom_positions();
            const double* p = hostPositions_.data();
            for (int i = 0; i < selectedCount_; ++i)
                host[i] = cvm::rvector(p[i], p[i+selectedCount_], p[i+2*selectedCount_]);
        }
        const double energy = proxy_->evaluateColvars();
        if (includeForces && selectedCount_ > 0) {
            if (!proxy_->has_gpu_support()) {
                const auto& host = *proxy_->get_atom_applied_forces();
                double* f = hostForces_.data();
                for (int i = 0; i < selectedCount_; ++i) {
                    f[i] = host[i].x;
                    f[i+selectedCount_] = host[i].y;
                    f[i+2*selectedCount_] = host[i].z;
                }
                checkCuda(cudaMemcpyAsync(reinterpret_cast<void*>(forces_.getDevicePointer()), f,
                                          3*selectedCount_*sizeof(double), cudaMemcpyHostToDevice, stream),
                          "uploading selected forces");
            }
            // Colvars group streams enqueue waits on this stream before returning.
            // Consequently scatter sees completed force sums without a host barrier.
            scatter_->execute(selectedCount_);
        }
        return includeEnergy ? energy : 0.0;
    }

    colvarproxy_openmm& proxy() override { return *proxy_; }
    std::string backend() const override { return backend_; }

private:
    void updateMapping() {
        if (selectedCount_ == 0) return;
        const auto& ids = proxy_->requestedAtoms();
        if (mapping_->dirty) {
            const auto& order = cuda_.getAtomIndex();
            for (int i = 0; i < cuda_.getNumAtoms(); ++i)
                inverseOrder_[order[i]] = i;
            for (int i = 0; i < selectedCount_; ++i)
                reordered_[i] = inverseOrder_[ids[i]];
            selected_.upload(reordered_);
        }
        // setPositions can clear cell offsets without causing a reorder.  Check
        // only M selected offsets on the host; upload only when they change.
        bool offsetsChanged = mapping_->dirty;
        const auto& offsets = cuda_.getPosCellOffsets();
        for (int i = 0; i < selectedCount_; ++i) {
            const auto& next = offsets[reordered_[i]];
            auto& old = cellOffsets_[i];
            if (old.x != next.x || old.y != next.y || old.z != next.z) {
                old = next;
                offsetsChanged = true;
            }
        }
        if (offsetsChanged) offsets_.upload(cellOffsets_);
        mapping_->dirty = false;
    }

    void describeBackend() {
        std::ostringstream out;
        out << "CUDA/" << (proxy_->has_gpu_support() ? "ColvarsGPU" : "ColvarsCPU")
            << "; selected atoms=" << selectedCount_
            << "; precision=" << (cuda_.getUseDoublePrecision() ? "double" :
                                  cuda_.getUseMixedPrecision() ? "mixed" : "single");
        if (proxy_->has_gpu_support()) {
            out << "; GPU atom groups; CPU bias updates; CPU CVCs=";
            bool first = true;
            for (const auto* variable : *proxy_->module().variables()) {
                for (const auto& component : variable->get_cvcs()) {
                    if (!component->is_enabled(colvardeps::f_cvc_support_gpu)) {
                        if (!first) out << ',';
                        out << variable->name << '/' << component->function_type();
                        first = false;
                    }
                }
            }
            if (first) out << "none";
        } else {
            out << "; per force evaluation D2H=" << 3*sizeof(double)*selectedCount_
                << " bytes positions, H2D=" << 3*sizeof(double)*selectedCount_ << " bytes forces";
        }
        backend_ = out.str();
    }

    CudaContext& cuda_;
    std::shared_ptr<MappingState> mapping_;
    std::unique_ptr<CudaProxy> proxy_;
    CudaArray positions_, forces_, selected_, offsets_;
    ComputeKernel gather_, scatter_;
    std::vector<int> reordered_, inverseOrder_;
    std::vector<mm_int4> cellOffsets_;
    PinnedBuffer hostPositions_, hostForces_;
    int selectedCount_ = 0;
    std::string backend_;
};

} // namespace

std::unique_ptr<CudaColvarsKernel> ColvarsPlugin::createCudaColvarsKernel(
    ContextImpl& context, const ColvarsForce& owner) {
    auto* data = reinterpret_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
    if (data == nullptr || data->contexts.size() != 1)
        throw OpenMMException("Colvars CUDA currently requires exactly one OpenMM CUDA device");
    return std::make_unique<ResidentCudaColvarsKernel>(context, owner, *data->contexts[0]);
}
