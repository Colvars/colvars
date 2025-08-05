#include "colvar_gpu_support.h"
#include "colvarmodule.h"

namespace colvars_gpu {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
int gpuAssert(cudaError_t code, const char *file, int line)
{
  if (code != cudaSuccess) {
    const std::string error =
      std::string("GPUassert: ") +
      cudaGetErrorString(code) + file;
    return cvm::error(error, COLVARS_ERROR);
  }
  return COLVARS_OK;
}
#endif
}
