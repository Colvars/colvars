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

class colvarproxy_gpu {
public:
  colvarproxy_gpu():
    support_gpu(false),
    gpu_code_type_used(gpu_code_t::CPU),
    gpu_id(0) {
    gpu_code_type_supported.push_back(gpu_code_t::CPU);
#if defined(COLVARS_CUDA)
    gpu_code_type_supported.push_back(gpu_code_t::CUDA);
    support_gpu = true;
#elif defined(COLVARS_HIP)
    gpu_code_type_supported.push_back(gpu_code_type::HIP);
    support_gpu = true;
#elif defined(COLVARS_SYCL)
    gpu_code_type_supported.push_back(gpu_code_type::SYCL);
    support_gpu = true;
#endif // COLVARS_SYCL
  }
  virtual int set_gpu(gpu_dev_id_t* gpu_id_in = nullptr);
  bool has_gpu_support() const {
    return support_gpu;
  }
  gpu_code_t get_gpu_code_type() const {
    return gpu_code_type_used;
  }
  int set_gpu_code_type(gpu_code_t gpu_code_type_in);
  virtual int create_stream(gpu_stream_t* stream, gpu_dev_id_t* gpu_id_in = nullptr);
  virtual int sync_all_streams();
  virtual int get_default_device(gpu_dev_id_t* device) const;
  template <typename T>
  int allocate_device(T **pp, const size_t len, gpu_dev_id_t* gpu_id_in = nullptr) {
    return allocate_device_T((void **)pp, len, sizeof(T), gpu_id_in);
  }
  template <typename T>
  int allocate_device_async(T **pp, const size_t len, gpu_stream_t* stream, gpu_dev_id_t* gpu_id_in = nullptr) {
    return allocate_device_T_async((void **)pp, len, sizeof(T), stream, gpu_id_in);
  }
  template <typename T>
  int deallocate_device(T **pp) {
    return deallocate_device_T((void **)pp);
  }
  template <typename T>
  int deallocate_device_async(T **pp, gpu_stream_t* stream) {
    return deallocate_device_T_async((void **)pp, stream);
  }
  virtual int allocate_device_T(void **pp, const size_t len, const size_t sizeofT, gpu_dev_id_t* gpu_id_in = nullptr);
  virtual int deallocate_device_T(void **pp);
  virtual int allocate_device_T_async(void **pp, const size_t len, const size_t sizeofT, gpu_stream_t* stream, gpu_dev_id_t* gpu_id_in = nullptr);
  virtual int deallocate_device_T_async(void **pp, gpu_stream_t* stream);
  ~colvarproxy_gpu();
protected:
  bool support_gpu;
  gpu_code_t gpu_code_type_used;
  std::vector<gpu_code_t> gpu_code_type_supported;
  std::vector<gpu_stream_t> gpu_streams;
  gpu_dev_id_t gpu_id;
};

} // namespace colvars_gpu

#endif // COLVAR_GPU_SUPPORT_H
