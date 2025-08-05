#include "colvarproxy_gpu.h"
#include "colvarmodule.h"

using namespace colvars_gpu;

int colvarproxy_gpu::allocate_host_T(void **pp, const size_t len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaMallocHost(pp, sizeofT*len));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::deallocate_host_T(void **pp) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  if (*pp != nullptr) {
    error_code |= checkGPUError(cudaFreeHost((void *)(*pp)));
    *pp = nullptr;
  }
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::allocate_device_T(void **pp, const size_t len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
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

int colvarproxy_gpu::allocate_device_T_async(void **pp, const size_t len, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaMallocAsync(pp, sizeofT*len, stream));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::deallocate_device_T_async(void **pp, cudaStream_t stream) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  if (*pp != nullptr) {
    error_code |= checkGPUError(cudaFreeAsync((void *)(*pp), stream));
    *pp = nullptr;
  }
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::clear_device_array_T(void *data, const size_t ndata, const size_t sizeofT) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  // if (data != nullptr) {
  error_code |= checkGPUError(cudaMemset(data, 0, sizeofT*ndata));
  // }
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::clear_device_array_T_async(void *data, const size_t ndata, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  // if (data != nullptr) {
  error_code |= checkGPUError(cudaMemsetAsync(data, 0, sizeofT*ndata, stream));
  // }
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::copy_HtoD_T(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaMemcpy(d_array, h_array, sizeofT*array_len, cudaMemcpyHostToDevice));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::copy_HtoD_T_async(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaMemcpyAsync(d_array, h_array, sizeofT*array_len, cudaMemcpyHostToDevice, stream));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::copy_DtoH_T(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaMemcpy(h_array, d_array, sizeofT*array_len, cudaMemcpyDeviceToHost));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::copy_DtoH_T_async(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaMemcpyAsync(h_array, d_array, sizeofT*array_len, cudaMemcpyDeviceToHost, stream));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::copy_DtoD_T(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaMemcpy(d_dst, d_src, sizeofT*array_len, cudaMemcpyDeviceToDevice));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

int colvarproxy_gpu::copy_DtoD_T_async(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  error_code |= checkGPUError(cudaMemcpyAsync(d_dst, d_src, sizeofT*array_len, cudaMemcpyDeviceToDevice, stream));
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
  error_code = COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}

colvarproxy_gpu::~colvarproxy_gpu() {
}
