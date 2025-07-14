#ifndef COLVARPROXY_GPU_H
#define COLVARPROXY_GPU_H

#include "colvar_gpu_support.h"
#include "colvarmodule.h"

class colvarproxy_gpu {
public:
  colvarproxy_gpu():
    support_gpu(false),
    gpu_code_type_used(colvars_gpu::gpu_code_t::CPU),
    gpu_id(0) {
    gpu_code_type_supported.push_back(colvars_gpu::gpu_code_t::CPU);
#if defined(COLVARS_CUDA)
    gpu_code_type_supported.push_back(colvars_gpu::gpu_code_t::CUDA);
    support_gpu = true;
#elif defined(COLVARS_HIP)
    gpu_code_type_supported.push_back(gpu_code_type::HIP);
    support_gpu = true;
#elif defined(COLVARS_SYCL)
    gpu_code_type_supported.push_back(gpu_code_type::SYCL);
    support_gpu = true;
#endif // COLVARS_SYCL
  }
  virtual int set_gpu(colvars_gpu::gpu_dev_id_t* gpu_id_in = nullptr);
  bool has_gpu_support() const {
    return support_gpu;
  }
  colvars_gpu::gpu_code_t get_gpu_code_type() const {
    return gpu_code_type_used;
  }
  int set_gpu_code_type(colvars_gpu::gpu_code_t gpu_code_type_in);
  virtual int create_stream(colvars_gpu::gpu_stream_t* stream, colvars_gpu::gpu_dev_id_t* gpu_id_in = nullptr);
  virtual int sync_all_streams();
  virtual int get_default_device(colvars_gpu::gpu_dev_id_t* device) const;
  template <typename T>
  int allocate_device(T **pp, const size_t len, colvars_gpu::gpu_dev_id_t* gpu_id_in = nullptr) {
    return allocate_device_T((void **)pp, len, sizeof(T), gpu_id_in);
  }
  template <typename T>
  int allocate_device_async(T **pp, const size_t len, colvars_gpu::gpu_stream_t* stream, colvars_gpu::gpu_dev_id_t* gpu_id_in = nullptr) {
    return allocate_device_T_async((void **)pp, len, sizeof(T), stream, gpu_id_in);
  }
  template <typename T>
  int deallocate_device(T **pp) {
    return deallocate_device_T((void **)pp);
  }
  template <typename T>
  int deallocate_device_async(T **pp, colvars_gpu::gpu_stream_t* stream) {
    return deallocate_device_T_async((void **)pp, stream);
  }
  virtual int allocate_device_T(void **pp, const size_t len, const size_t sizeofT, colvars_gpu::gpu_dev_id_t* gpu_id_in = nullptr);
  virtual int deallocate_device_T(void **pp);
  virtual int allocate_device_T_async(void **pp, const size_t len, const size_t sizeofT, colvars_gpu::gpu_stream_t* stream, colvars_gpu::gpu_dev_id_t* gpu_id_in = nullptr);
  virtual int deallocate_device_T_async(void **pp, colvars_gpu::gpu_stream_t* stream);
  virtual cvm::real* proxy_atoms_masses_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_charges_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_positions_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_total_forces_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_new_colvar_forces_gpu() {return nullptr;}
  ~colvarproxy_gpu();
protected:
  bool support_gpu;
  colvars_gpu::gpu_code_t gpu_code_type_used;
  std::vector<colvars_gpu::gpu_code_t> gpu_code_type_supported;
  std::vector<colvars_gpu::gpu_stream_t> gpu_streams;
  colvars_gpu::gpu_dev_id_t gpu_id;
};

#endif // COLVARPROXY_GPU_H
