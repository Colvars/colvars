#include "colvar_gpu_support.h"
#include "colvarmodule.h"

#if defined(__has_include)
# if __has_include(<stacktrace>)
#  define use_cpp_stacktrace
#  include <stacktrace>
# endif

#endif

namespace colvars_gpu {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
int gpuAssert(cudaError_t code, const char *file, int line)
{
  if (code != cudaSuccess) {
    std::string error =
      std::string("GPUassert: ") +
      cudaGetErrorString(code) + " " + file + ", line " + cvm::to_str(line);
#ifdef use_cpp_stacktrace
#ifdef __cpp_lib_stacktrace
      error += "\nBacktrace:\n" + std::to_string(std::stacktrace::current()) + "\n";
#endif
#endif
    return cvm::error(error, COLVARS_ERROR);
  }
  return COLVARS_OK;
}

int add_clear_array_node_impl(
  void* dst, const size_t num_elements, const size_t sizeofT,
  cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
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
  if (cvm::debug()) {
    cvm::log(
      "Add a memset clear node: ptr = " + cvm::to_str(dst) +
      " width = " + cvm::to_str(width) + " elementSize = " +
      cvm::to_str(elementSize) + "\n");
  }
  return checkGPUError(cudaGraphAddMemsetNode(
    &node_out, graph, dependencies.data(),
    dependencies.size(), &memsetParams));
}

int add_copy_node_impl(
  const void* src, void* dst, const size_t num_elements, const size_t sizeofT,
  cudaMemcpyKind kind, cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  // cudaMemcpy3DParms    memcpyParams     = {0};
  // memcpyParams.kind     = kind;
  // memcpyParams.srcArray = NULL;
  // memcpyParams.srcPos   = make_cudaPos(0, 0, 0);
  // memcpyParams.srcPtr   = make_cudaPitchedPtr(
  //   (void*)src, sizeofT * num_elements, num_elements, 1);
  // memcpyParams.dstArray = NULL;
  // memcpyParams.dstPos   = make_cudaPos(0, 0, 0);
  // memcpyParams.dstPtr   = make_cudaPitchedPtr(
  //   (void*)dst, sizeofT * num_elements, num_elements, 1);
  // memcpyParams.extent   = make_cudaExtent(sizeofT * num_elements, 1, 1);
  if (cvm::debug()) {
    cvm::log(
      "Add a memcpy node: src = " + cvm::to_str(src) +
      " dst = " + cvm::to_str(dst) + " num_elements = " +
      cvm::to_str(num_elements) + "\n");
  }
  // return checkGPUError(
  //   cudaGraphAddMemcpyNode(
  //     &node_out, graph, dependencies.data(),
  //     dependencies.size(), &memcpyParams));
  return checkGPUError(cudaGraphAddMemcpyNode1D(
    &node_out, graph, dependencies.data(), dependencies.size(),
    dst, src, sizeofT * num_elements, kind));
}
#endif
}
