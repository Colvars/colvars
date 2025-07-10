#include "colvar_gpu_support.h"
#include "colvarmodule.h"

namespace colvars_gpu {

int gpuAssert(gpu_error_t code, const char *file, int line)
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

int colvarproxy_gpu::set_gpu(gpu_dev_id_t* gpu_id_in) {
  if (gpu_id_in == nullptr) {
    return get_default_device(&gpu_id);
  } else {
    gpu_id = *gpu_id_in;
    return COLVARS_OK;
  }
}

int colvarproxy_gpu::set_gpu_code_type(gpu_code_t gpu_code_type_in) {
  for (auto it = gpu_code_type_supported.begin();
         it != gpu_code_type_supported.end(); ++it) {
    if ((*it) == gpu_code_type_in) {
      gpu_code_type_used = gpu_code_type_in;
      return COLVARS_OK;
    }
  }
  return cvm::error("Unsupported GPU code type.", COLVARS_NOT_IMPLEMENTED);
}

int colvarproxy_gpu::create_stream(gpu_stream_t* stream, gpu_dev_id_t* gpu_id_in) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  if (gpu_id_in != nullptr) {
    if ((*gpu_id_in) != gpu_id) {
      error_code |= checkGPUError(cudaSetDevice(*gpu_id_in));
    }
  }
  error_code |= checkGPUError(cudaStreamCreate(stream));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::sync_all_streams() {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  for (auto it = gpu_streams.begin(); it != gpu_streams.end(); ++it) {
    error_code |= checkGPUError(cudaStreamSynchronize(*it));
  }
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::get_default_device(gpu_dev_id_t* device) const {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaGetDevice(device));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::allocate_device_T(void **pp, const size_t len, const size_t sizeofT, gpu_dev_id_t* gpu_id_in) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  if (gpu_id_in != nullptr) {
    if ((*gpu_id_in) != gpu_id) {
      error_code |= checkGPUError(cudaSetDevice(*gpu_id_in));
    }
  }
  error_code |= checkGPUError(cudaMalloc(pp, sizeofT*len));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::deallocate_device_T(void **pp) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  if (*pp != nullptr) {
    error_code |= checkGPUError(cudaFree((void *)(*pp)));
    *pp = nullptr;
  }
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::allocate_device_T_async(void **pp, const size_t len, const size_t sizeofT, gpu_stream_t* stream, gpu_dev_id_t* gpu_id_in) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  if (gpu_id_in != nullptr) {
    if ((*gpu_id_in) != gpu_id) {
      error_code |= checkGPUError(cudaSetDevice(*gpu_id_in));
    }
  }
  error_code |= checkGPUError(cudaMallocAsync(pp, sizeofT*len, *stream));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::deallocate_device_T_async(void **pp, gpu_stream_t* stream) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  if (*pp != nullptr) {
    error_code |= checkGPUError(cudaFreeAsync((void *)(*pp), *stream));
    *pp = nullptr;
  }
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

colvarproxy_gpu::~colvarproxy_gpu() {
  sync_all_streams();
  for (auto it = gpu_streams.begin(); it != gpu_streams.end(); ++it) {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
    checkGPUError(cudaStreamDestroy((*it)));
#elif defined(COLVARS_SYCL)
    // TODO: SYCL
#endif
  }
}

}
