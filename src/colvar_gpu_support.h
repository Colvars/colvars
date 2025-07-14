#ifndef COLVAR_GPU_SUPPORT_H
#define COLVAR_GPU_SUPPORT_H

#include <vector>

#if defined(COLVARS_CUDA)
#include <cuda_runtime.h>
#endif // defined(COLVARS_CUDA)

#if defined(COLVARS_HIP)
#include <hip/hip_runtime.h>

#ifndef cudaError_t
#define cudaError_t hipError_t
#endif

#ifndef cudaSuccess
#define cudaSuccess hipSuccess
#endif

#ifndef cudaGetErrorString
#define cudaGetErrorString hipGetErrorString
#endif

#ifndef cudaStreamDestroy
#define cudaStreamDestroy hipStreamDestroy
#endif

#ifndef cudaStreamCreate
#define cudaStreamCreate hipStreamCreate
#endif

#ifndef cudaStreamSynchronize
#define cudaStreamSynchronize hipStreamSynchronize
#endif

#ifndef cudaStream_t
#define cudaStream_t hipStream_t
#endif

#ifndef cudaSetDevice
#define cudaSetDevice hipSetDevice
#endif

#ifndef cudaGetDevice
#define cudaGetDevice hipGetDevice
#endif

#ifndef cudaMalloc
#define cudaMalloc hipMalloc
#endif

#ifndef cudaFree
#define cudaFree hipFree
#endif

#ifndef cudaMallocAsync
#define cudaMallocAsync hipMallocAsync
#endif

#ifndef cudaFreeAsync
#define cudaFreeAsync hipFreeAsync
#endif

#ifndef cudaHostAllocMapped
#define cudaHostAllocMapped hipHostMallocMapped
#endif

#ifndef cudaHostAlloc
#define cudaHostAlloc hipHostMalloc
#endif

#ifndef cudaFreeHost
#define cudaFreeHost hipHostFree
#endif
#endif // defined(COLVARS_HIP)

namespace colvars_gpu {

enum class gpu_code_t {
  CUDA, HIP, SYCL, CPU
};

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
#define COLVARS_HOST_DEVICE __device__ __host__
#else
#define COLVARS_HOST_DEVICE
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


#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
using gpu_stream_t = cudaStream_t;
using gpu_error_t = cudaError_t;
using gpu_dev_id_t = int;
#elif defined(COLVARS_SYCL)
using gpu_stream_t = sycl::queue;
#else
using gpu_stream_t = int;
using gpu_error_t = int;
using gpu_dev_id_t = int;
#endif

#if defined(COLVARS_CUDA) || defined (COLVARS_HIP)
#define checkGPUError(ans) gpuAssert((ans), __FILE__, __LINE__);
int gpuAssert(gpu_error_t code, const char *file, int line);
#endif

} // namespace colvars_gpu

#endif // COLVAR_GPU_SUPPORT_H
