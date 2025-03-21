#ifndef COLVARPROXY_CUDAGLOBALMASTER_KERNEL_H
#define COLVARPROXY_CUDAGLOBALMASTER_KERNEL_H

#include "colvartypes.h"
#include <cuda_runtime.h>

/**
 * @brief: Convert the device data from xxxyyyzzz to xyzxyzxyz and copy them to the host memory
 */
void transpose_to_host_rvector(
  const double* d_data_in,
  cvm::rvector* h_data_out,
  const int num_atoms,
  cudaStream_t stream);

/**
 * @brief: Convert the device data from xyzxyzxyz to xxxyyyzzz and copy them back to the device memory
 */
void transpose_from_host_rvector(
  double* d_data_out,
  const cvm::rvector* h_data_in,
  const int num_atoms,
  cudaStream_t stream);

/**
 * @brief: Convert the device data from float to double and copy them to the host memory
 */
void copy_float_to_host_double(
  const float* d_data_in,
  cvm::real* h_data_out,
  const int num_atoms,
  cudaStream_t stream);

#endif // COLVARPROXY_CUDAGLOBALMASTER_KERNEL_H
