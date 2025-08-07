#ifndef COLVARPROXY_GPU_H
#define COLVARPROXY_GPU_H

#include "colvar_gpu_support.h"
#include "colvarmodule.h"

class colvarproxy_gpu {
public:
  colvarproxy_gpu(): support_gpu(false) {}
  bool has_gpu_support() const {
    return support_gpu;
  }
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP) || defined (COLVARS_SYCL)
  virtual cudaStream_t get_default_stream() {return (cudaStream_t)0;}
  template <typename T>
  int allocate_host(T **pp, const size_t len) {
    return allocate_host_T((void **)pp, len, sizeof(T));
  }
  template <typename T>
  int deallocate_host(T **pp) {
    return deallocate_host_T((void **)pp);
  }
  template <typename T>
  int allocate_device(T **pp, const size_t len) {
    return allocate_device_T((void **)pp, len, sizeof(T));
  }
  template <typename T>
  int reallocate_device(T **pp, const size_t len) {
    int error_code = COLVARS_OK;
    error_code |= deallocate_device(pp);
    error_code |= allocate_device_T((void **)pp, len, sizeof(T));
    return error_code;
  }
  template <typename T>
  int reallocate_host(T **pp, const size_t len) {
    int error_code = COLVARS_OK;
    error_code |= deallocate_host(pp);
    error_code |= allocate_host_T((void **)pp, len, sizeof(T));
    return error_code;
  }
  template <typename T>
  int allocate_device_async(T **pp, const size_t len, cudaStream_t stream) {
    return allocate_device_T_async((void **)pp, len, sizeof(T), stream);
  }
  template <typename T>
  int deallocate_device(T **pp) {
    return deallocate_device_T((void **)pp);
  }
  template <typename T>
  int deallocate_device_async(T **pp, cudaStream_t stream) {
    return deallocate_device_T_async((void **)pp, stream);
  }
  template <typename T>
  int clear_device_array(T *data, const size_t ndata) {
    return clear_device_array_T(data, ndata, sizeof(T));
  }
  template <typename T>
  int clear_device_array_async(T *data, const size_t ndata, cudaStream_t stream) {
    return clear_device_array_T_async(data, ndata, sizeof(T), stream);
  }
  template <typename T>
  int copy_HtoD(const T *h_array, T *d_array, size_t array_len) {
    return copy_HtoD_T(h_array, d_array, array_len, sizeof(T));
  }
  template <typename T>
  int copy_HtoD_async(const T *h_array, T *d_array, size_t array_len, cudaStream_t stream) {
    return copy_HtoD_T_async(h_array, d_array, array_len, sizeof(T), stream);
  }
  template <typename T>
  int copy_DtoH(const T *d_array, T *h_array, size_t array_len) {
    return copy_DtoH_T(d_array, h_array, array_len, sizeof(T));
  }
  template <typename T>
  int copy_DtoH_async(const T *d_array, T *h_array, size_t array_len, cudaStream_t stream) {
    return copy_DtoH_T_async(d_array, h_array, array_len, sizeof(T), stream);
  }
  template <typename T>
  int copy_DtoD(const T *d_src, T *d_dst, size_t array_len) {
    return copy_DtoD_T(d_src, d_dst, array_len, sizeof(T));
  }
  template <typename T>
  int copy_DtoD_async(const T *d_src, T *d_dst, size_t array_len, cudaStream_t stream) {
    return copy_DtoD_T_async(d_src, d_dst, array_len, sizeof(T), stream);
  }
  virtual int allocate_host_T(void **pp, const size_t len, const size_t sizeofT);
  virtual int deallocate_host_T(void **pp);
  virtual int allocate_device_T(void **pp, const size_t len, const size_t sizeofT);
  virtual int deallocate_device_T(void **pp);
  virtual int clear_device_array_T(void *data, const size_t ndata, const size_t sizeofT);
  virtual int allocate_device_T_async(void **pp, const size_t len, const size_t sizeofT, cudaStream_t stream);
  virtual int deallocate_device_T_async(void **pp, cudaStream_t stream);
  virtual int clear_device_array_T_async(void *data, const size_t ndata, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_HtoD_T(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT);
  virtual int copy_HtoD_T_async(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_DtoH_T(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT);
  virtual int copy_DtoH_T_async(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_DtoD_T(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT);
  virtual int copy_DtoD_T_async(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  virtual float* proxy_atoms_masses_gpu_float() {return nullptr;}
  virtual float* proxy_atoms_charges_gpu_float() {return nullptr;}
  virtual cvm::real* proxy_atoms_masses_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_charges_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_positions_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_total_forces_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_new_colvar_forces_gpu() {return nullptr;}
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  virtual ~colvarproxy_gpu();
protected:
  bool support_gpu;
};

#endif // COLVARPROXY_GPU_H
