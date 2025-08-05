#include "colvar_gpu_support.h"
#include "colvaratoms_kernel.h"
#include "colvartypes.h"

#if defined(COLVARS_CUDA)
#include <cub/block/block_reduce.cuh>
#endif

// TODO: HIP CUB

namespace colvars_gpu {
#if defined(COLVARS_CUDA) || defined(COVLARS_HIP)
__global__ void atoms_pos_from_proxy_kernel(
  const int* __restrict atoms_proxy_index,
  const cvm::real* __restrict atoms_pos_x_proxy,
  const cvm::real* __restrict atoms_pos_y_proxy,
  const cvm::real* __restrict atoms_pos_z_proxy,
  cvm::real* __restrict atoms_pos_x_ag,
  cvm::real* __restrict atoms_pos_y_ag,
  cvm::real* __restrict atoms_pos_z_ag,
  int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    const int proxy_index = atoms_proxy_index[i];
    atoms_pos_x_ag[i] = atoms_pos_x_proxy[proxy_index];
    atoms_pos_y_ag[i] = atoms_pos_y_proxy[proxy_index];
    atoms_pos_z_ag[i] = atoms_pos_z_proxy[proxy_index];
  }
}

int atoms_pos_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_pos_proxy,
  cvm::real* atoms_pos_ag,
  int num_atoms,
  int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_pos_x_proxy = atoms_pos_proxy;
  const cvm::real* atoms_pos_y_proxy = atoms_pos_x_proxy + proxy_stride;
  const cvm::real* atoms_pos_z_proxy = atoms_pos_y_proxy + proxy_stride;
  cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] =
    {&atoms_proxy_index,
     &atoms_pos_x_proxy,
     &atoms_pos_y_proxy,
     &atoms_pos_z_proxy,
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)atoms_pos_from_proxy_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

int atoms_pos_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_pos_proxy,
  cvm::real* atoms_pos_ag,
  int num_atoms,
  int proxy_stride,
  cudaStream_t stream) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_pos_x_proxy = atoms_pos_proxy;
  const cvm::real* atoms_pos_y_proxy = atoms_pos_x_proxy + proxy_stride;
  const cvm::real* atoms_pos_z_proxy = atoms_pos_y_proxy + proxy_stride;
  cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] =
    {&atoms_proxy_index,
     &atoms_pos_x_proxy,
     &atoms_pos_y_proxy,
     &atoms_pos_z_proxy,
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &num_atoms};
  return checkGPUError(cudaLaunchKernel((void*)atoms_pos_from_proxy_kernel,
    num_blocks, block_size, args, 0, stream));
}

__global__ void change_one_coordinate_kernel(
  cvm::real* __restrict atoms_pos_ag,
  size_t array_id, cvm::real step_size) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i == 0) {
    atoms_pos_ag[array_id] += step_size;
  }
}

int change_one_coordinate(
  cvm::real* atoms_pos_ag, size_t atom_id_in_group, int xyz,
  cvm::real step_size, size_t num_atoms, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  if (xyz > 0 && xyz < 3) {
    size_t array_id = num_atoms * xyz + atom_id_in_group;
    void* args[] = {&atoms_pos_ag, &array_id, &step_size};
    error_code |= checkGPUError(cudaLaunchKernel(
      (void*)change_one_coordinate_kernel,
      1, 1, args, 0, stream));
  }
  return error_code;
}

int change_one_coordinate(
  cvm::real* atoms_pos_ag,
  size_t atom_id_in_group, int xyz,
  cvm::real step_size,
  size_t num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  size_t array_id = num_atoms * xyz + atom_id_in_group;
  void* args[] = {&atoms_pos_ag, &array_id, &step_size};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)change_one_coordinate_kernel;
  kernelNodeParams.gridDim        = dim3(1, 1, 1);
  kernelNodeParams.blockDim       = dim3(1, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (node == 0) {
    return checkGPUError(cudaGraphAddKernelNode(
      &node, graph, dependencies.data(),
      dependencies.size(), &kernelNodeParams));
  } else {
    return checkGPUError(cudaGraphKernelNodeSetParams(
      node, &kernelNodeParams));
  }
}

template <int BLOCK_SIZE, bool b_cog, bool b_com>
__global__ void atoms_calc_cog_com_kernel(
  const cvm::real* __restrict atoms_mass,
  const cvm::real* __restrict atoms_pos_x_ag,
  const cvm::real* __restrict atoms_pos_y_ag,
  const cvm::real* __restrict atoms_pos_z_ag,
  cvm::rvector* __restrict cog_out,
  cvm::rvector* __restrict com_out,
  cvm::rvector* __restrict h_cog_out,
  cvm::rvector* __restrict h_com_out,
  cvm::real total_mass,
  unsigned int* __restrict tbcount,
  int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ bool isLastBlockDone;
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  cvm::rvector cog{0, 0, 0};
  cvm::rvector com{0, 0, 0};
  if (i < num_atoms) {
    if (b_cog) {
      cog.x = atoms_pos_x_ag[i];
      cog.y = atoms_pos_y_ag[i];
      cog.z = atoms_pos_z_ag[i];
    }
    if (b_com) {
      com.x = atoms_mass[i] * atoms_pos_x_ag[i];
      com.y = atoms_mass[i] * atoms_pos_y_ag[i];
      com.z = atoms_mass[i] * atoms_pos_z_ag[i];
    }
  }
  __syncthreads();
  typedef cub::BlockReduce<double, BLOCK_SIZE> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  if (b_cog) {
    cog.x = BlockReduce(temp_storage).Sum(cog.x); __syncthreads();
    cog.y = BlockReduce(temp_storage).Sum(cog.y); __syncthreads();
    cog.z = BlockReduce(temp_storage).Sum(cog.z); __syncthreads();
  }
  if (b_com) {
    com.x = BlockReduce(temp_storage).Sum(com.x); __syncthreads();
    com.y = BlockReduce(temp_storage).Sum(com.y); __syncthreads();
    com.z = BlockReduce(temp_storage).Sum(com.z); __syncthreads();
  }
  if (threadIdx.x == 0) {
    if (b_cog) {
      atomicAdd(&(cog_out->x), cog.x);
      atomicAdd(&(cog_out->y), cog.y);
      atomicAdd(&(cog_out->z), cog.z);
    }
    if (b_com) {
      atomicAdd(&(com_out->x), com.x);
      atomicAdd(&(com_out->y), com.y);
      atomicAdd(&(com_out->z), com.z);
    }
    __threadfence();
    unsigned int value = atomicInc(tbcount, gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      if (b_cog) {
        h_cog_out->x = cog_out->x;
        h_cog_out->y = cog_out->y;
        h_cog_out->z = cog_out->z;
      }
      if (b_com) {
        com.x = com_out->x / total_mass;
        com.y = com_out->y / total_mass;
        com.z = com_out->z / total_mass;
        com_out->x = com.x;
        com_out->y = com.y;
        com_out->z = com.z;
        h_com_out->x = com.x;
        h_com_out->y = com.y;
        h_com_out->z = com.z;
      }
    }
  }
}

int atoms_calc_cog_com(
  const cvm::real* atoms_pos_ag,
  const cvm::real* atoms_mass,
  int num_atoms,
  cvm::rvector* cog_out,
  cvm::rvector* com_out,
  cvm::rvector* h_cog_out,
  cvm::rvector* h_com_out,
  cvm::real total_mass,
  unsigned int* tbcount,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  const cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  const cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] =
    {&atoms_mass,
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &cog_out,
     &com_out,
     &h_cog_out,
     &h_com_out,
     &total_mass,
     &tbcount,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)atoms_calc_cog_com_kernel<block_size, true, true>;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

int atoms_calc_cog(
  const cvm::real* atoms_pos_ag,
  int num_atoms,
  cvm::rvector* cog_out,
  cvm::rvector* h_cog_out,
  unsigned int* tbcount,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  const cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  const cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  const cvm::real* atoms_mass = nullptr;
  cvm::real total_mass = 0.0;
  void* args[] =
    {&atoms_mass,
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &cog_out,
     &h_cog_out,
     &total_mass,
     &tbcount,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)atoms_calc_cog_com_kernel<block_size, true, false>;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

template <bool ag_rotate>
__global__ void atoms_total_force_from_proxy_kernel(
  const int* __restrict atoms_proxy_index,
  const cvm::real* __restrict atoms_total_force_x_proxy,
  const cvm::real* __restrict atoms_total_force_y_proxy,
  const cvm::real* __restrict atoms_total_force_z_proxy,
  cvm::real* __restrict atoms_total_force_x_ag,
  cvm::real* __restrict atoms_total_force_y_ag,
  cvm::real* __restrict atoms_total_force_z_ag,
  cvm::quaternion* __restrict q,
  int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  cvm::rmatrix rot_mat;
  if (ag_rotate) {
    rot_mat = q->rotation_matrix();
  }
  if (i < num_atoms) {
    const int proxy_index = atoms_proxy_index[i];
    const cvm::real fx = atoms_total_force_x_proxy[proxy_index];
    const cvm::real fy = atoms_total_force_y_proxy[proxy_index];
    const cvm::real fz = atoms_total_force_z_proxy[proxy_index];
    if (ag_rotate) {
      atoms_total_force_x_ag[i] = rot_mat.xx * fx + rot_mat.xy * fy + rot_mat.xz * fz;
      atoms_total_force_y_ag[i] = rot_mat.yx * fx + rot_mat.yy * fy + rot_mat.yz * fz;
      atoms_total_force_z_ag[i] = rot_mat.zx * fx + rot_mat.zy * fy + rot_mat.zz * fz;
    } else {
      atoms_total_force_x_ag[i] = fx;
      atoms_total_force_y_ag[i] = fy;
      atoms_total_force_z_ag[i] = fz;
    }
  }
}

void atoms_total_force_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_total_force_proxy,
  cvm::real* atoms_total_force_ag,
  bool rotate,
  cvm::quaternion* q,
  int num_atoms,
  int proxy_stride,
  cudaStream_t stream) {
  if (num_atoms == 0) return;
  const int block_size = default_block_size;
  const int grid = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_total_force_x_proxy = atoms_total_force_proxy;
  const cvm::real* atoms_total_force_y_proxy = atoms_total_force_x_proxy + proxy_stride;
  const cvm::real* atoms_total_force_z_proxy = atoms_total_force_y_proxy + proxy_stride;
  cvm::real* atoms_total_force_x_ag = atoms_total_force_ag;
  cvm::real* atoms_total_force_y_ag = atoms_total_force_x_ag + num_atoms;
  cvm::real* atoms_total_force_z_ag = atoms_total_force_y_ag + num_atoms;
  void* args[] =
    {&atoms_proxy_index,
     &atoms_total_force_x_proxy,
     &atoms_total_force_y_proxy,
     &atoms_total_force_z_proxy,
     &atoms_total_force_x_ag,
     &atoms_total_force_y_ag,
     &atoms_total_force_z_ag,
     &q,
     &num_atoms};
  if (rotate) {
    checkGPUError(cudaLaunchKernel(
      (void*)atoms_total_force_from_proxy_kernel<true>,
      grid, block_size, args, 0, stream));
  } else {
    checkGPUError(cudaLaunchKernel(
      (void*)atoms_total_force_from_proxy_kernel<false>,
      grid, block_size, args, 0, stream));
  }
}

template <bool ag_rotate>
__global__ void apply_colvar_force_to_proxy_kernel(
  const int* __restrict atoms_proxy_index,
  cvm::real* __restrict atoms_applied_force_x_proxy,
  cvm::real* __restrict atoms_applied_force_y_proxy,
  cvm::real* __restrict atoms_applied_force_z_proxy,
  const cvm::real* __restrict grad_x,
  const cvm::real* __restrict grad_y,
  const cvm::real* __restrict grad_z,
  cvm::real* force_ptr,
  cvm::quaternion* q,
  int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  const cvm::real force = (*force_ptr);
  const cvm::rmatrix rot_inv = q->conjugate().rotation_matrix();
  if (i < num_atoms) {
    const int proxy_index = atoms_proxy_index[i];
    cvm::real fx, fy, fz;
    if (ag_rotate) {
      fx = force * (rot_inv.xx * grad_x[i] +
                    rot_inv.xy * grad_y[i] +
                    rot_inv.xz * grad_z[i]);
      fy = force * (rot_inv.yx * grad_x[i] +
                    rot_inv.yy * grad_y[i] +
                    rot_inv.yz * grad_z[i]);
      fz = force * (rot_inv.zx * grad_x[i] +
                    rot_inv.zy * grad_y[i] +
                    rot_inv.zz * grad_z[i]);
    } else {
      fx = force * grad_x[i];
      fy = force * grad_y[i];
      fz = force * grad_z[i];
    }
    atomicAdd(&(atoms_applied_force_x_proxy[proxy_index]), fx);
    atomicAdd(&(atoms_applied_force_y_proxy[proxy_index]), fy);
    atomicAdd(&(atoms_applied_force_z_proxy[proxy_index]), fz);
  }
}

int apply_colvar_force_to_proxy(
  const int* atoms_proxy_index,
  cvm::real* atoms_applied_force_proxy,
  const cvm::real* atoms_grad_ag,
  cvm::real* colvar_force,
  bool rotate,
  cvm::quaternion* q,
  int num_atoms,
  int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  // if (num_atoms == 0) return;
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  cvm::real* atoms_applied_force_x_proxy = atoms_applied_force_proxy;
  cvm::real* atoms_applied_force_y_proxy = atoms_applied_force_x_proxy + proxy_stride;
  cvm::real* atoms_applied_force_z_proxy = atoms_applied_force_y_proxy + proxy_stride;
  const cvm::real* grad_x = atoms_grad_ag;
  const cvm::real* grad_y = grad_x + num_atoms;
  const cvm::real* grad_z = grad_y + num_atoms;
  void* args[] = {
    &atoms_proxy_index,
    &atoms_applied_force_x_proxy,
    &atoms_applied_force_y_proxy,
    &atoms_applied_force_z_proxy,
    &grad_x,
    &grad_y,
    &grad_z,
    &colvar_force,
    &q,
    &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (rotate) {
    // checkGPUError(cudaLaunchKernel(
    //   (void*)apply_colvar_force_to_proxy_kernel<true>,
    //   num_blocks, block_size, args, 0, stream));
    kernelNodeParams.func = (void*)apply_colvar_force_to_proxy_kernel<true>;
  } else {
    // checkGPUError(cudaLaunchKernel(
    //   (void*)apply_colvar_force_to_proxy_kernel<false>,
    //   num_blocks, block_size, args, 0, stream));
    kernelNodeParams.func = (void*)apply_colvar_force_to_proxy_kernel<false>;
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

/**
 * @brief Calculate the gradients of rotated positions wrt quaternion
 */
template <bool B_ag_center,
          bool B_ag_rotate,
          int BLOCK_SIZE>
__global__ void calc_fit_forces_impl_loop1_kernel(
  const cvm::real* __restrict atoms_grad_or_force_x,
  const cvm::real* __restrict atoms_grad_or_force_y,
  const cvm::real* __restrict atoms_grad_or_force_z,
  const cvm::real* __restrict atoms_pos_unrotated_x,
  const cvm::real* __restrict atoms_pos_unrotated_y,
  const cvm::real* __restrict atoms_pos_unrotated_z,
  const cvm::quaternion* __restrict q,
  const int main_group_size,
  const int group_for_fit_size,
  double3* __restrict atom_grad,
  double4* __restrict sum_dxdq,
  unsigned int* __restrict count) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ bool isLastBlockDone;
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  double sum_dxdq_x = 0;
  double sum_dxdq_y = 0;
  double sum_dxdq_z = 0;
  double sum_dxdq_w = 0;
  cvm::rvector main_grad{0, 0, 0};
  if (i < main_group_size) {
    if (B_ag_center || B_ag_rotate) {
      main_grad = cvm::rvector{
        atoms_grad_or_force_x[i],
        atoms_grad_or_force_y[i],
        atoms_grad_or_force_z[i]};
    }
    if (B_ag_rotate) {
      cvm::quaternion const dxdq =
        q->position_derivative_inner(
          cvm::rvector{
            atoms_pos_unrotated_x[i],
            atoms_pos_unrotated_y[i],
            atoms_pos_unrotated_z[i]},
          main_grad);
      sum_dxdq_x = dxdq[0];
      sum_dxdq_y = dxdq[1];
      sum_dxdq_z = dxdq[2];
      sum_dxdq_w = dxdq[3];
    }
  }
  __syncthreads();
  typedef cub::BlockReduce<double, BLOCK_SIZE> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  main_grad.x = BlockReduce(temp_storage).Sum(main_grad.x); __syncthreads();
  main_grad.y = BlockReduce(temp_storage).Sum(main_grad.y); __syncthreads();
  main_grad.z = BlockReduce(temp_storage).Sum(main_grad.z); __syncthreads();
  sum_dxdq_x = BlockReduce(temp_storage).Sum(sum_dxdq_x); __syncthreads();
  sum_dxdq_y = BlockReduce(temp_storage).Sum(sum_dxdq_y); __syncthreads();
  sum_dxdq_z = BlockReduce(temp_storage).Sum(sum_dxdq_z); __syncthreads();
  sum_dxdq_w = BlockReduce(temp_storage).Sum(sum_dxdq_w); __syncthreads();
  if (threadIdx.x == 0) {
    atomicAdd(&(atom_grad->x), main_grad.x);
    atomicAdd(&(atom_grad->y), main_grad.y);
    atomicAdd(&(atom_grad->z), main_grad.z);
    atomicAdd(&(sum_dxdq->x), sum_dxdq_x);
    atomicAdd(&(sum_dxdq->y), sum_dxdq_y);
    atomicAdd(&(sum_dxdq->z), sum_dxdq_z);
    atomicAdd(&(sum_dxdq->w), sum_dxdq_w);
    __threadfence();
    unsigned int value = atomicInc(count, gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      if (B_ag_center) {
        main_grad.x = atom_grad->x;
        main_grad.y = atom_grad->y;
        main_grad.z = atom_grad->z;
        if (B_ag_rotate) {
          const cvm::rmatrix rot_inv = q->conjugate().rotation_matrix();
          const cvm::real x = main_grad.x * rot_inv.xx + main_grad.y * rot_inv.xy + main_grad.z * rot_inv.xz;
          const cvm::real y = main_grad.x * rot_inv.yx + main_grad.y * rot_inv.yy + main_grad.z * rot_inv.yz;
          const cvm::real z = main_grad.x * rot_inv.zx + main_grad.y * rot_inv.zy + main_grad.z * rot_inv.zz;
          main_grad.x = x;
          main_grad.y = y;
          main_grad.z = z;
        }
        main_grad.x *= -1.0 / group_for_fit_size;
        main_grad.y *= -1.0 / group_for_fit_size;
        main_grad.z *= -1.0 / group_for_fit_size;
        atom_grad->x = main_grad.x;
        atom_grad->y = main_grad.y;
        atom_grad->z = main_grad.z;
      }
    }
  }
}

int calc_fit_gradients_impl_loop1(
  const cvm::real* pos_unrotated,
  cvm::real* main_grad,
  const cvm::quaternion* q,
  int num_atoms_main,
  int num_atoms_fitting,
  double3* atom_grad,
  double4* sum_dxdq,
  unsigned int* tbcount,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms_main + block_size - 1) / block_size;
  const cvm::real* pos_x = pos_unrotated;
  const cvm::real* pos_y = pos_x + num_atoms_main;
  const cvm::real* pos_z = pos_y + num_atoms_main;
  const cvm::real* grad_x = main_grad;
  const cvm::real* grad_y = grad_x + num_atoms_main;
  const cvm::real* grad_z = grad_y + num_atoms_main;
  // auto access_fitting = [grad_x, grad_y, grad_z] __device__ (int i){
  //   return cvm::rvector{grad_x[i], grad_y[i], grad_z[i]};
  // };
  void* args[] = {
    // &access_fitting,
    &grad_x, &grad_y, &grad_z,
    &pos_x, &pos_y, &pos_z,
    &q, &num_atoms_main, &num_atoms_fitting,
    &atom_grad, &sum_dxdq, &tbcount};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop1_kernel<true, true, block_size>;
  }
  if (ag_center && !ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop1_kernel<true, false, block_size>;
  }
  if (!ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop1_kernel<false, true, block_size>;
  }
  if (!ag_center && !ag_rotate) {
    return COLVARS_OK;
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

// loop 2: iterate over the fitting group
template <bool B_ag_center,
          bool B_ag_rotate,
          int BLOCK_SIZE,
          bool do_fit_gradients>
__global__ void calc_fit_forces_impl_loop2_kernel(
  cvm::real* __restrict fit_grad_x,
  cvm::real* __restrict fit_grad_y,
  cvm::real* __restrict fit_grad_z,
  colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const double3* __restrict atom_grad,
  const double4* __restrict sum_dxdq,
  const int* __restrict atoms_proxy_index,
  cvm::real* __restrict proxy_new_force_x,
  cvm::real* __restrict proxy_new_force_y,
  cvm::real* __restrict proxy_new_force_z,
  const int group_for_fit_size) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < group_for_fit_size) {
    cvm::rvector fitting_force_grad{0, 0, 0};
    if (B_ag_center) {
      fitting_force_grad.x += atom_grad->x;
      fitting_force_grad.y += atom_grad->y;
      fitting_force_grad.z += atom_grad->z;
    }
    if (B_ag_rotate) {
      cvm::rvector dq0_1[4];
      #pragma unroll
      for (int j = 0; j < 4; ++j) {dq0_1[j].set(0);}
      rot_deriv->calc_derivative_wrt_group1<false, true>(i, nullptr, dq0_1);
      fitting_force_grad += sum_dxdq->x * dq0_1[0] +
                            sum_dxdq->y * dq0_1[1] +
                            sum_dxdq->z * dq0_1[2] +
                            sum_dxdq->w * dq0_1[3];
    }
    // accessor_fitting(i, fitting_force_grad);
    if (do_fit_gradients) {
      fit_grad_x[i] = fitting_force_grad.x;
      fit_grad_y[i] = fitting_force_grad.y;
      fit_grad_z[i] = fitting_force_grad.z;
    } else {
      const int pid = atoms_proxy_index[i];
      atomicAdd(&(proxy_new_force_x[pid]), fitting_force_grad.x);
      atomicAdd(&(proxy_new_force_y[pid]), fitting_force_grad.y);
      atomicAdd(&(proxy_new_force_z[pid]), fitting_force_grad.z);
    }
  }
}

int calc_fit_gradients_impl_loop2(
  cvm::real* fit_grad,
  colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const double3* atom_grad,
  const double4* sum_dxdq,
  int group_for_fit_size,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (group_for_fit_size + block_size - 1) / block_size;
  cvm::real* fit_grad_x = fit_grad;
  cvm::real* fit_grad_y = fit_grad_x + group_for_fit_size;
  cvm::real* fit_grad_z = fit_grad_y + group_for_fit_size;
  // auto access_fitting = [fit_grad_x, fit_grad_y, fit_grad_z] __device__ (
  //   int i, const cvm::rvector& grad){
  //   fit_grad_x[i] = grad.x;
  //   fit_grad_y[i] = grad.y;
  //   fit_grad_z[i] = grad.z;
  // };
  const int* atoms_proxy_index = nullptr;
  cvm::real* proxy_new_force_x = nullptr;
  cvm::real* proxy_new_force_y = nullptr;
  cvm::real* proxy_new_force_z = nullptr;
  void* args[] = {
    // &access_fitting,
    &fit_grad_x, &fit_grad_y, &fit_grad_z,
    &rot_deriv, &atom_grad, &sum_dxdq,
    &atoms_proxy_index,
    &proxy_new_force_x,
    &proxy_new_force_y,
    &proxy_new_force_z,
    &group_for_fit_size};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<true, true, block_size, true>;
  }
  if (ag_center && !ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<true, false, block_size, true>;
  }
  if (!ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<false, true, block_size, true>;
  }
  if (!ag_center && !ag_rotate) {
    return COLVARS_OK;
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

__global__ void apply_translation_kernel(
  cvm::real* __restrict atoms_pos_x_ag,
  cvm::real* __restrict atoms_pos_y_ag,
  cvm::real* __restrict atoms_pos_z_ag,
  cvm::real translation_vector_factor,
  const cvm::rvector* __restrict translation_vector,
  int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    atoms_pos_x_ag[i] += translation_vector_factor * translation_vector->x;
    atoms_pos_y_ag[i] += translation_vector_factor * translation_vector->y;
    atoms_pos_z_ag[i] += translation_vector_factor * translation_vector->z;
  }
}

int apply_translation(
  cvm::real* atoms_pos_ag,
  cvm::real translation_vector_factor,
  const cvm::rvector* translation_vector,
  int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] = {
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &translation_vector_factor,
     &translation_vector,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)apply_translation_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

__global__ void rotate_with_quaternion_kernel(
  cvm::real* __restrict pos_x,
  cvm::real* __restrict pos_y,
  cvm::real* __restrict pos_z,
  cvm::quaternion* __restrict q, int num_atoms) {
  const auto rot_mat = q->rotation_matrix();
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    const cvm::real new_x = rot_mat.xx * pos_x[i] +
                            rot_mat.xy * pos_y[i] +
                            rot_mat.xz * pos_z[i];
    const cvm::real new_y = rot_mat.yx * pos_x[i] +
                            rot_mat.yy * pos_y[i] +
                            rot_mat.yz * pos_z[i];
    const cvm::real new_z = rot_mat.zx * pos_x[i] +
                            rot_mat.zy * pos_y[i] +
                            rot_mat.zz * pos_z[i];
    pos_x[i] = new_x;
    pos_y[i] = new_y;
    pos_z[i] = new_z;
  }
}

int rotate_with_quaternion(
  cvm::real* atoms_pos_ag,
  cvm::quaternion* q,
  int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] = {
    &atoms_pos_x_ag,
    &atoms_pos_y_ag,
    &atoms_pos_z_ag,
    &q, &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)rotate_with_quaternion_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

__global__ void accumulate_cpu_force_kernel(
  const cvm::real* __restrict h_atoms_force_x,
  const cvm::real* __restrict h_atoms_force_y,
  const cvm::real* __restrict h_atoms_force_z,
  cvm::real* __restrict d_atoms_force_x,
  cvm::real* __restrict d_atoms_force_y,
  cvm::real* __restrict d_atoms_force_z,
  const int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    d_atoms_force_x[i] += h_atoms_force_x[i];
    d_atoms_force_y[i] += h_atoms_force_y[i];
    d_atoms_force_z[i] += h_atoms_force_z[i];
  }
}

int accumulate_cpu_force(
  const cvm::real* h_atoms_force,
  cvm::real* d_atoms_force,
  int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const cvm::real* h_atoms_force_x = h_atoms_force;
  const cvm::real* h_atoms_force_y = h_atoms_force_x + num_atoms;
  const cvm::real* h_atoms_force_z = h_atoms_force_y + num_atoms;
  cvm::real* d_atoms_force_x = d_atoms_force;
  cvm::real* d_atoms_force_y = d_atoms_force_x + num_atoms;
  cvm::real* d_atoms_force_z = d_atoms_force_y + num_atoms;
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  void* args[] = {
    &h_atoms_force_x,
    &h_atoms_force_y,
    &h_atoms_force_z,
    &d_atoms_force_x,
    &d_atoms_force_y,
    &d_atoms_force_z,
    &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)accumulate_cpu_force;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

__global__ void apply_force_with_inverse_rotation_kernel(
  const cvm::real* __restrict atoms_force_x,
  const cvm::real* __restrict atoms_force_y,
  const cvm::real* __restrict atoms_force_z,
  const cvm::quaternion* q,
  int* __restrict atoms_proxy_index,
  cvm::real* __restrict proxy_new_force_x,
  cvm::real* __restrict proxy_new_force_y,
  cvm::real* __restrict proxy_new_force_z,
  const int num_atoms) {
  const cvm::rmatrix rot_inv = q->conjugate().rotation_matrix();
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    const cvm::rvector f_ia{
      rot_inv.xx * atoms_force_x[i] +
      rot_inv.xy * atoms_force_y[i] +
      rot_inv.xz * atoms_force_z[i],
      rot_inv.yx * atoms_force_x[i] +
      rot_inv.yy * atoms_force_y[i] +
      rot_inv.yz * atoms_force_z[i],
      rot_inv.zx * atoms_force_x[i] +
      rot_inv.zy * atoms_force_y[i] +
      rot_inv.zz * atoms_force_z[i],
    };
    const int pid = atoms_proxy_index[i];
    atomicAdd(&(proxy_new_force_x[pid]), f_ia.x);
    atomicAdd(&(proxy_new_force_y[pid]), f_ia.y);
    atomicAdd(&(proxy_new_force_z[pid]), f_ia.z);
  }
}

int apply_force_with_inverse_rotation(
  const cvm::real* atoms_force,
  const cvm::quaternion* q,
  const int* atoms_proxy_index,
  cvm::real* proxy_new_force,
  int num_atoms,
  int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_force_x = atoms_force;
  const cvm::real* atoms_force_y = atoms_force_x + num_atoms;
  const cvm::real* atoms_force_z = atoms_force_y + num_atoms;
  cvm::real* proxy_new_force_x = proxy_new_force;
  cvm::real* proxy_new_force_y = proxy_new_force_x + proxy_stride;
  cvm::real* proxy_new_force_z = proxy_new_force_y + proxy_stride;
  void* args[] =
    {&atoms_force_x,
     &atoms_force_y,
     &atoms_force_z,
     &q,
     &atoms_proxy_index,
     &proxy_new_force_x,
     &proxy_new_force_y,
     &proxy_new_force_z,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)apply_force_with_inverse_rotation_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

__global__ void apply_force_kernel(
  const cvm::real* __restrict atoms_force_x,
  const cvm::real* __restrict atoms_force_y,
  const cvm::real* __restrict atoms_force_z,
  int* __restrict atoms_proxy_index,
  cvm::real* __restrict proxy_new_force_x,
  cvm::real* __restrict proxy_new_force_y,
  cvm::real* __restrict proxy_new_force_z,
  const int num_atoms) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    const int pid = atoms_proxy_index[i];
    atomicAdd(&(proxy_new_force_x[pid]), atoms_force_x[i]);
    atomicAdd(&(proxy_new_force_y[pid]), atoms_force_y[i]);
    atomicAdd(&(proxy_new_force_z[pid]), atoms_force_z[i]);
  }
}

int apply_force(
  const cvm::real* atoms_force,
  const int* atoms_proxy_index,
  cvm::real* proxy_new_force,
  int num_atoms,
  int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_force_x = atoms_force;
  const cvm::real* atoms_force_y = atoms_force_x + num_atoms;
  const cvm::real* atoms_force_z = atoms_force_y + num_atoms;
  cvm::real* proxy_new_force_x = proxy_new_force;
  cvm::real* proxy_new_force_y = proxy_new_force_x + proxy_stride;
  cvm::real* proxy_new_force_z = proxy_new_force_y + proxy_stride;
  void* args[] =
    {&atoms_force_x,
     &atoms_force_y,
     &atoms_force_z,
     &atoms_proxy_index,
     &proxy_new_force_x,
     &proxy_new_force_y,
     &proxy_new_force_z,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)apply_force_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

int calc_fit_forces_impl_loop1(
  const cvm::real* pos_unrotated,
  cvm::real* main_force,
  const cvm::quaternion* q,
  int num_atoms_main,
  int num_atoms_fitting,
  double3* atom_grad,
  double4* sum_dxdq,
  unsigned int* tbcount,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  return calc_fit_gradients_impl_loop1(
    pos_unrotated, main_force, q, num_atoms_main,
    num_atoms_fitting, atom_grad, sum_dxdq, tbcount,
    ag_center, ag_rotate, node, graph, dependencies);
}

int calc_fit_forces_impl_loop2(
  colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const double3* atom_grad,
  const double4* sum_dxdq,
  const int* atoms_proxy_index,
  cvm::real* proxy_new_force,
  int group_for_fit_size,
  int proxy_stride,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (group_for_fit_size + block_size - 1) / block_size;
  cvm::real* fit_grad_x = nullptr;
  cvm::real* fit_grad_y = nullptr;
  cvm::real* fit_grad_z = nullptr;
  cvm::real* proxy_new_force_x = proxy_new_force;
  cvm::real* proxy_new_force_y = proxy_new_force_x + proxy_stride;
  cvm::real* proxy_new_force_z = proxy_new_force_y + proxy_stride;
  void* args[] = {
    // &access_fitting,
    &fit_grad_x, &fit_grad_y, &fit_grad_z,
    &rot_deriv, &atom_grad, &sum_dxdq,
    &atoms_proxy_index,
    &proxy_new_force_x,
    &proxy_new_force_y,
    &proxy_new_force_z,
    &group_for_fit_size};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<true, true, block_size, false>;
  }
  if (ag_center && !ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<true, false, block_size, false>;
  }
  if (!ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<false, true, block_size, false>;
  }
  if (!ag_center && !ag_rotate) {
    return COLVARS_OK;
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

#elif defined(COLVARS_SYCL)
#endif // defined(COLVARS_CUDA) || defined(COVLARS_HIP)
};
