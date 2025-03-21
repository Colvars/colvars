#include "colvarproxy_cudaglobalmaster_kernel.h"

__global__ void transpose_to_host_rvector_kernel(
  const double* __restrict d_data_in,
  cvm::rvector* __restrict h_data_out,
  const int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    // printf("i = %d\n", i);
    h_data_out[i].x = d_data_in[i];
    h_data_out[i].y = d_data_in[i + num_atoms];
    h_data_out[i].z = d_data_in[i + num_atoms * 2];
    // printf("i = %d, x = %lf, y = %lf, z = %lf\n", i,
    //        d_data_in[i], d_data_in[i + num_atoms], d_data_in[i + num_atoms * 2]);
  }
}

void transpose_to_host_rvector(
  const double* d_data_in,
  cvm::rvector* h_data_out,
  const int num_atoms,
  cudaStream_t stream) {
  const int block_size = 128;
  const int grid = (num_atoms + block_size - 1) / block_size;
  if (grid == 0) return;
  transpose_to_host_rvector_kernel<<<grid, block_size, 0, stream>>>(
    d_data_in, h_data_out, num_atoms);
}

__global__ void transpose_from_host_rvector_kernel(
  double*             __restrict d_data_out,
  const cvm::rvector* __restrict h_data_in,
  const int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    d_data_out[i]                 = h_data_in[i].x;
    d_data_out[i + num_atoms]     = h_data_in[i].y;
    d_data_out[i + num_atoms * 2] = h_data_in[i].z;
  }
}

void transpose_from_host_rvector(
  double* d_data_out,
  const cvm::rvector* h_data_in,
  const int num_atoms,
  cudaStream_t stream) {
  const int block_size = 128;
  const int grid = (num_atoms + block_size - 1) / block_size;
  if (grid == 0) return;
  transpose_from_host_rvector_kernel<<<grid, block_size, 0, stream>>>(
    d_data_out, h_data_in, num_atoms);
}

__global__ void copy_float_to_host_double_kernel(
  const float*  __restrict d_data_in,
  cvm::real*    __restrict h_data_out,
  const int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    h_data_out[i] = cvm::real(d_data_in[i]);
  }
}

void copy_float_to_host_double(
  const float* d_data_in,
  cvm::real* h_data_out,
  const int num_atoms,
  cudaStream_t stream) {
  const int block_size = 128;
  const int grid = (num_atoms + block_size - 1) / block_size;
  if (grid == 0) return;
  copy_float_to_host_double_kernel<<<grid, block_size, 0, stream>>>(
    d_data_in, h_data_out, num_atoms);
}
