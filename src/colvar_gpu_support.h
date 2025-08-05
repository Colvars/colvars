#ifndef COLVAR_GPU_SUPPORT_H
#define COLVAR_GPU_SUPPORT_H

#include <vector>

#define COLVARS_STRINGIFY(s) STRINGIFY_HELPER(s)
#define STRINGIFY_HELPER(s) #s

#if defined(COLVARS_CUDA)
#include <cuda_runtime.h>
#endif // defined(COLVARS_CUDA)

#if defined(COLVARS_HIP)
#ifndef cudaError_t
#define cudaError_t hipError_t
#endif // cudaError_t

#ifndef cudaFree
#define cudaFree hipFree
#endif // cudaFree

#ifndef cudaFreeHost
#define cudaFreeHost hipFreeHost
#endif // cudaFreeHost

#ifndef cudaGetErrorString
#define cudaGetErrorString hipGetErrorString
#endif // cudaGetErrorString

#ifndef cudaGraphAddChildGraphNode
#define cudaGraphAddChildGraphNode hipGraphAddChildGraphNode
#endif // cudaGraphAddChildGraphNode

#ifndef cudaGraphAddKernelNode
#define cudaGraphAddKernelNode hipGraphAddKernelNode
#endif // cudaGraphAddKernelNode

#ifndef cudaGraphAddMemcpyNode
#define cudaGraphAddMemcpyNode hipGraphAddMemcpyNode
#endif // cudaGraphAddMemcpyNode

#ifndef cudaGraphAddMemsetNode
#define cudaGraphAddMemsetNode hipGraphAddMemsetNode
#endif // cudaGraphAddMemsetNode

#ifndef cudaGraphCreate
#define cudaGraphCreate hipGraphCreate
#endif // cudaGraphCreate

#ifndef cudaGraphDestroy
#define cudaGraphDestroy hipGraphDestroy
#endif // cudaGraphDestroy

#ifndef cudaGraphExecDestroy
#define cudaGraphExecDestroy hipGraphExecDestroy
#endif // cudaGraphExecDestroy

#ifndef cudaGraphExecMemcpyNodeSetParams
#define cudaGraphExecMemcpyNodeSetParams hipGraphExecMemcpyNodeSetParams
#endif // cudaGraphExecMemcpyNodeSetParams

#ifndef cudaGraphExec_t
#define cudaGraphExec_t hipGraphExec_t
#endif // cudaGraphExec_t

#ifndef cudaGraphInstantiate
#define cudaGraphInstantiate hipGraphInstantiate
#endif // cudaGraphInstantiate

#ifndef cudaGraphLaunch
#define cudaGraphLaunch hipGraphLaunch
#endif // cudaGraphLaunch

#ifndef cudaGraphNode_t
#define cudaGraphNode_t hipGraphNode_t
#endif // cudaGraphNode_t

#ifndef cudaGraph_t
#define cudaGraph_t hipGraph_t
#endif // cudaGraph_t

#ifndef cudaHostAllocMapped
#define cudaHostAllocMapped hipHostAllocMapped
#endif // cudaHostAllocMapped

#ifndef cudaKernelNodeParams
#define cudaKernelNodeParams hipKernelNodeParams
#endif // cudaKernelNodeParams

#ifndef cudaMalloc
#define cudaMalloc hipMalloc
#endif // cudaMalloc

#ifndef cudaMallocHost
#define cudaMallocHost hipMallocHost
#endif // cudaMallocHost

#ifndef cudaMemcpy
#define cudaMemcpy hipMemcpy
#endif // cudaMemcpy

#ifndef cudaMemcpy3DParms
#define cudaMemcpy3DParms hipMemcpy3DParms
#endif // cudaMemcpy3DParms

#ifndef cudaMemcpyAsync
#define cudaMemcpyAsync hipMemcpyAsync
#endif // cudaMemcpyAsync

#ifndef cudaMemcpyDeviceToDevice
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#endif // cudaMemcpyDeviceToDevice

#ifndef cudaMemcpyDeviceToHost
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#endif // cudaMemcpyDeviceToHost

#ifndef cudaMemcpyHostToDevice
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#endif // cudaMemcpyHostToDevice

#ifndef cudaMemcpyKind
#define cudaMemcpyKind hipMemcpyKind
#endif // cudaMemcpyKind

#ifndef cudaMemset
#define cudaMemset hipMemset
#endif // cudaMemset

#ifndef cudaMemsetAsync
#define cudaMemsetAsync hipMemsetAsync
#endif // cudaMemsetAsync

#ifndef cudaStreamCreate
#define cudaStreamCreate hipStreamCreate
#endif // cudaStreamCreate

#ifndef cudaStreamDestroy
#define cudaStreamDestroy hipStreamDestroy
#endif // cudaStreamDestroy

#ifndef cudaStreamSynchronize
#define cudaStreamSynchronize hipStreamSynchronize
#endif // cudaStreamSynchronize

#ifndef cudaStream_t
#define cudaStream_t hipStream_t
#endif // cudaStream_t

#ifndef cudaSuccess
#define cudaSuccess hipSuccess
#endif // cudaSuccess

#ifndef make_cudaExtent
#define make_cudaExtent make_hipExtent
#endif // make_cudaExtent

#ifndef make_cudaPitchedPtr
#define make_cudaPitchedPtr make_hipPitchedPtr
#endif // make_cudaPitchedPtr

#ifndef make_cudaPos
#define make_cudaPos make_hipPos
#endif // make_cudaPos

#endif // defined(COLVARS_HIP)

namespace colvars_gpu {

constexpr int default_block_size = 128;

// enum class gpu_code_t {
//   CUDA, HIP, SYCL, CPU
// };

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
#define COLVARS_HOST_DEVICE __device__ __host__
#define COLVARS_DEVICE __device__
#else
#define COLVARS_HOST_DEVICE
#define COLVARS_DEVICE
#endif

// TODO: What about SYCL?
#if ( defined(COLVARS_CUDA) || defined(COLVARS_HIP) )
template <typename T>
class CudaHostAllocator {
public:
  using value_type = T;

  CudaHostAllocator() = default;

  template<typename U>
  constexpr CudaHostAllocator(const CudaHostAllocator<U>&) noexcept {}

  friend bool operator==(const CudaHostAllocator&, const CudaHostAllocator&) { return true; }
  friend bool operator!=(const CudaHostAllocator&, const CudaHostAllocator&) { return false; }

  T* allocate(size_t n) {
    T* ptr;
    if (cudaHostAlloc(&ptr, n * sizeof(T), cudaHostAllocMapped) != cudaSuccess) {
      throw std::bad_alloc();
    }
    return ptr;
  }
  void deallocate(T* ptr, size_t n) noexcept {
    cudaFreeHost(ptr);
  }
  template<typename U, typename... Args>
  void construct(U* p, Args&&... args) {
      new(p) U(std::forward<Args>(args)...);
  }

  template<typename U>
  void destroy(U* p) noexcept {
      p->~U();
  }
};
#endif


#if defined(COLVARS_CUDA) || defined (COLVARS_HIP)
int gpuAssert(cudaError_t code, const char *file, int line);
#endif

} // namespace colvars_gpu

#if defined(COLVARS_CUDA) || defined (COLVARS_HIP)
#define checkGPUError(ans) colvars_gpu::gpuAssert((ans), __FILE__, __LINE__);
#endif

namespace colvars_gpu {
#if defined(COLVARS_CUDA) || defined (COLVARS_HIP)
template <typename T>
int add_clear_array_node(
  T* dst, const size_t num_elements,
  cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  // size_t elementSize, width;
  const size_t sizeofT = sizeof(T);
  /**< Size of each element in bytes. Must be 1, 2, or 4. */
  const size_t elementSize =
    (sizeofT % 4 == 0) ? 4 :
    ((sizeofT % 2 == 0) ? 2 : 1);
  const size_t width = num_elements * (sizeofT / elementSize);
  cudaMemsetParams memsetParams = {0};
  memsetParams.dst         = (void*)dst;
  memsetParams.value       = 0;
  memsetParams.elementSize = elementSize;
  memsetParams.width       = width;
  memsetParams.height      = 1;
  return checkGPUError(cudaGraphAddMemsetNode(
    &node_out, graph, dependencies.data(),
    dependencies.size(), &memsetParams));
}

template <typename T>
int add_copy_node(
  const T* src, T* dst, size_t num_elements,
  cudaMemcpyKind kind, cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  cudaMemcpy3DParms    memcpyParams     = {0};
  memcpyParams.kind     = kind;
  memcpyParams.srcArray = NULL;
  memcpyParams.srcPos   = make_cudaPos(0, 0, 0);
  memcpyParams.srcPtr   = make_cudaPitchedPtr(
    (void*)src, sizeof(T) * num_elements, num_elements, 1);
  memcpyParams.dstArray = NULL;
  memcpyParams.dstPos   = make_cudaPos(0, 0, 0);
  memcpyParams.dstPtr   = make_cudaPitchedPtr(
    (void*)dst, sizeof(T) * num_elements, num_elements, 1);
  memcpyParams.extent   = make_cudaExtent(sizeof(T) * num_elements, 1, 1);
  return checkGPUError(
    cudaGraphAddMemcpyNode(
      &node_out, graph, dependencies.data(),
      dependencies.size(), &memcpyParams));
}

template <typename T>
int update_copy_node(
  const T* src, T* dst, size_t num_elements,
  cudaMemcpyKind kind, cudaGraphNode_t& node,
  cudaGraphExec_t& graph_exec) {
  cudaMemcpy3DParms    memcpyParams     = {0};
  memcpyParams.kind     = kind;
  memcpyParams.srcArray = NULL;
  memcpyParams.srcPos   = make_cudaPos(0, 0, 0);
  memcpyParams.srcPtr   = make_cudaPitchedPtr(
    (void*)src, sizeof(T) * num_elements, num_elements, 1);
  memcpyParams.dstArray = NULL;
  memcpyParams.dstPos   = make_cudaPos(0, 0, 0);
  memcpyParams.dstPtr   = make_cudaPitchedPtr(
    (void*)dst, sizeof(T) * num_elements, num_elements, 1);
  memcpyParams.extent   = make_cudaExtent(sizeof(T) * num_elements, 1, 1);
  int error_code = checkGPUError(
    cudaGraphMemcpyNodeSetParams(node, &memcpyParams));
  error_code |= checkGPUError(
    cudaGraphExecMemcpyNodeSetParams(
    graph_exec, node,
    &memcpyParams));
  return error_code;
}
#endif
}

#define ADD_DEPENDENCY(fieldName, dependencies_vector, nodes_map) do {\
  const std::string s = COLVARS_STRINGIFY(fieldName) ;\
  try { dependencies_vector.push_back(nodes_map.at(s)); }\
  catch (const std::out_of_range& oor) { \
    return cvm::error(cvm::to_str("BUG: cannot find node ") + s); } \
} while (0);

#define ADD_DEPENDENCY_IF(fieldName, dependencies_vector, nodes_map) do {\
  const std::string s = COLVARS_STRINGIFY(fieldName) ;\
  if (auto search = nodes_map.find(s); search != nodes_map.end()) {\
    dependencies_vector.push_back(search->second);\
  }\
} while (0);

#endif // COLVAR_GPU_SUPPORT_H
