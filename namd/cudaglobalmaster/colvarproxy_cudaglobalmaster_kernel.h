#ifndef COLVARPROXY_CUDAGLOBALMASTER_KERNEL_H
#define COLVARPROXY_CUDAGLOBALMASTER_KERNEL_H

#include "colvartypes.h"

#if defined (NAMD_CUDA) || defined (NAMD_HIP)

#ifdef NAMD_CUDA
#include <cuda_runtime.h>
#endif  // NAMD_CUDA

#ifdef NAMD_HIP
#include <hip/hip_runtime.h>
#endif  // NAMD_HIP

#include "HipDefines.h"

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

#endif // defined (NAMD_CUDA) || defined (NAMD_HIP)

#endif // COLVARPROXY_CUDAGLOBALMASTER_KERNEL_H
