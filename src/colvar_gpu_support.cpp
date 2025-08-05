#include "colvar_gpu_support.h"
#include "colvarmodule.h"

namespace colvars_gpu {

int gpuAssert(cudaError_t code, const char *file, int line)
{
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  if (code != cudaSuccess) {
    const std::string error =
      std::string("GPUassert: ") +
      cudaGetErrorString(code) + file;
    return cvm::error(error, COLVARS_ERROR);
  }
#endif
  return COLVARS_OK;
}

}
