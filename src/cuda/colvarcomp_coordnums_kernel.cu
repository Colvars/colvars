#include "colvarcomp_coordnums_kernel.h"
#include "colvarcomp_coordnums.h"

#if defined(COLVARS_CUDA)
#include <cub/block/block_reduce.cuh>
#include <cooperative_groups.h>
#endif

#if defined (COLVARS_HIP)
#include <hipcub/block/block_reduce.hpp>
#include <hip/hip_cooperative_groups.h>
#define cub hipcub
#endif

#if defined(COLVARS_CUDA) || defined (COLVARS_HIP)
namespace colvars_gpu {
template <int N, int M, int blockSize, int group2BatchSize,
          int numGroup2BatchesPerBlock, int flags>
__global__ void computeCoordinationNumberTwoGroupsCUDAKernel1(
  const cvm::real* __restrict pos1x,
  const cvm::real* __restrict pos1y,
  const cvm::real* __restrict pos1z,
  const cvm::real* __restrict pos2x,
  const cvm::real* __restrict pos2y,
  const cvm::real* __restrict pos2z,
  const unsigned int numAtoms1,
  const unsigned int numAtoms2,
  const int en, const int ed,
  const cvm::rvector inv_r0_vec,
  const cvm::rvector inv_r0sq_vec,
#if 0
  // TODO: Wait for https://github.com/Colvars/colvars/pull/919
  const cvm::system_boundary_conditions bc,
#endif
  cvm::real* __restrict gx1,
  cvm::real* __restrict gy1,
  cvm::real* __restrict gz1,
  cvm::real* __restrict gx2,
  cvm::real* __restrict gy2,
  cvm::real* __restrict gz2,
  const cvm::real pairlist_tol,
  const cvm::real pairlist_tol_l2_max,
  bool* __restrict pairlist,
  unsigned int* __restrict tbcount,
  cvm::real* __restrict coordnum_tmp,
  cvm::real* __restrict coordnum_out) {
#if 0
  // TODO: Optimize for static en and ed. Wait for https://github.com/Colvars/colvars/pull/926.
  constexpr const bool static_exponents = (N > 0) && (M > 0);
#endif
  constexpr const bool use_pairlist = flags & colvar::coordnum::ef_use_pairlist;
  constexpr const bool rebuild_pairlist = flags & colvar::coordnum::ef_rebuild_pairlist;
  constexpr const bool gradients = flags & colvar::coordnum::ef_gradients;
  constexpr const bool use_internal_pbc = flags & colvar::coordnum::ef_use_internal_pbc;
  static_assert(use_internal_pbc == true, "The CUDA kernel requires internal PBC.");
  static_assert(blockSize == group2BatchSize * numGroup2BatchesPerBlock, "blockSize != group2BatchSize * numGroup2BatchesPerBlock");
  // Shared memory buffers for atoms in group2
  __shared__ double3 shPosition[group2BatchSize];
  __shared__ double3 shJGrad[numGroup2BatchesPerBlock][group2BatchSize];
  extern __shared__ bool shPairlist_buffer[];
  bool (&shPairlist)[numGroup2BatchesPerBlock][group2BatchSize][blockSize] =
    *reinterpret_cast<bool (*)[numGroup2BatchesPerBlock][group2BatchSize][blockSize]>(shPairlist_buffer);
  __shared__ bool shJMask[group2BatchSize];
  __shared__ bool isLastBlockDone;
  // Total coordnum
  cvm::real ei = 0;
  // Number of blocks required to iterate over group1
  const unsigned int numBlocksInGroup1 = (numAtoms1 + blockSize - 1) / blockSize;
  // Number of blocks required to iterate over group2
  const unsigned int numBatchesInGroup2 = (numAtoms2 + group2BatchSize - 1) / group2BatchSize;
  const unsigned int group2WorkSize = numBatchesInGroup2 * group2BatchSize;
  const unsigned int group2BatchID = threadIdx.x / group2BatchSize;
  const unsigned int group2LaneID = threadIdx.x % group2BatchSize;
  for (unsigned int i = blockIdx.x; i < numBlocksInGroup1; i += gridDim.x) {
    const unsigned int tid = i * blockDim.x + threadIdx.x;
    // Load the atom i from group1
    const bool mask_i = tid < numAtoms1;
    const cvm::real x1 = mask_i ? pos1x[tid] : 0;
    const cvm::real y1 = mask_i ? pos1y[tid] : 0;
    const cvm::real z1 = mask_i ? pos1z[tid] : 0;
    double3 iGrad{0, 0, 0};
    // Load atom j from group2
    for (unsigned int k = 0; k < group2WorkSize; k += group2BatchSize) {
      const unsigned int j = k + group2LaneID;
      const bool mask_j = j < numAtoms2;
      if (group2BatchID == 0) {
        if (mask_j) {
          shPosition[group2LaneID].x = pos2x[j];
          shPosition[group2LaneID].y = pos2y[j];
          shPosition[group2LaneID].z = pos2z[j];
        }
        shJMask[group2LaneID] = mask_j;
      }
      if constexpr (gradients) {
        shJGrad[group2BatchID][group2LaneID].x = 0;
        shJGrad[group2BatchID][group2LaneID].y = 0;
        shJGrad[group2BatchID][group2LaneID].z = 0;
      }
      if (use_pairlist && !(rebuild_pairlist)) {
        #pragma unroll
        for (unsigned int t = 0; t < group2BatchSize; ++t) {
          const int jid = k + t;
          const bool mask_jid = jid < numAtoms2;
          shPairlist[group2BatchID][t][threadIdx.x] =
            (mask_i && mask_jid) ? pairlist[tid+jid*numAtoms1] : false;
        }
      }
      __syncthreads();
      for (unsigned int t = 0; t < group2BatchSize; ++t) {
        const unsigned int jid = t ^ group2LaneID;
        const bool mask_jid = shJMask[jid];
        if (mask_i && mask_jid) {
          const cvm::real x2 = shPosition[jid].x;
          const cvm::real y2 = shPosition[jid].y;
          const cvm::real z2 = shPosition[jid].z;
#if 0
          // TODO: Wait for https://github.com/Colvars/colvars/pull/919
          if constexpr (!(use_pairlist)) {
            ei += compute_pair_coordnum<flags>(
              inv_r0_vec, inv_r0sq_vec, en, ed,
              x1, y1, z1, x2, y2, z2,
              iGrad.x, iGrad.y, iGrad.z,
              shJGrad[group2BatchID][jid].x,
              shJGrad[group2BatchID][jid].y,
              shJGrad[group2BatchID][jid].z,
              pairlist_tol, pairlist_tol_l2_max, bc);
          } else {
            if constexpr (!(rebuild_pairlist)) {
              const bool within = shPairlist[group2BatchID][jid][threadIdx.x];
              if (within) {
                ei += compute_pair_coordnum<flags>(
                  inv_r0_vec, inv_r0sq_vec, en, ed,
                  x1, y1, z1, x2, y2, z2,
                  iGrad.x, iGrad.y, iGrad.z,
                  shJGrad[group2BatchID][jid].x,
                  shJGrad[group2BatchID][jid].y,
                  shJGrad[group2BatchID][jid].z,
                  pairlist_tol, pairlist_tol_l2_max, bc);
              }
            } else {
              const double f = compute_pair_coordnum<flags>(
                  inv_r0_vec, inv_r0sq_vec, en, ed,
                  x1, y1, z1, x2, y2, z2,
                  iGrad.x, iGrad.y, iGrad.z,
                  shJGrad[group2BatchID][jid].x,
                  shJGrad[group2BatchID][jid].y,
                  shJGrad[group2BatchID][jid].z,
                  pairlist_tol, pairlist_tol_l2_max, bc);
              shPairlist[group2BatchID][jid][threadIdx.x] = f > 0.0;
              ei += f;
            }
          }
#endif
        }
        __syncthreads();
      }
      if constexpr (use_pairlist && rebuild_pairlist) {
        #pragma unroll
        for (unsigned int t = 0; t < group2BatchSize; ++t) {
          const int jid = k + t;
          const bool mask_jid = jid < numAtoms2;
          if (mask_i && mask_jid) {
            pairlist[tid+jid*numAtoms1] = shPairlist[group2BatchID][t][threadIdx.x];
          }
        }
      }
      if constexpr (gradients) {
        // Reduction over the shared memory
        #pragma unroll
        for (unsigned int l = numGroup2BatchesPerBlock / 2; l > 0; l >>= 1) {
          if (group2BatchID < l) {
            shJGrad[group2BatchID][group2LaneID].x += shJGrad[group2BatchID + l][group2LaneID].x;
            shJGrad[group2BatchID][group2LaneID].y += shJGrad[group2BatchID + l][group2LaneID].y;
            shJGrad[group2BatchID][group2LaneID].z += shJGrad[group2BatchID + l][group2LaneID].z;
          }
          __syncthreads();
        }
        if (group2BatchID == 0) {
          if (shJMask[group2LaneID]) {
            atomicAdd(&gx2[j], shJGrad[0][group2LaneID].x);
            atomicAdd(&gy2[j], shJGrad[0][group2LaneID].y);
            atomicAdd(&gz2[j], shJGrad[0][group2LaneID].z);
          }
        }
      }
      __syncthreads();
    }
    if constexpr (gradients) {
      if (mask_i) {
        // Save the i-gradients to group1
        atomicAdd(&gx1[tid], iGrad.x);
        atomicAdd(&gy1[tid], iGrad.y);
        atomicAdd(&gz1[tid], iGrad.z);
      }
    }
  }
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  // Reduction for coordnum
  typedef cub::BlockReduce<cvm::real, blockSize> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  const cvm::real total_e = BlockReduce(temp_storage).Sum(ei); __syncthreads();
  if (threadIdx.x == 0) {
    atomicAdd(coordnum_tmp, total_e);
    __threadfence();
    unsigned int value = atomicInc(tbcount, gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      *coordnum_out = *coordnum_tmp;
      *coordnum_tmp = 0;
      tbcount[0] = 0;
    }
  }
}

int calc_value_coordnum_two_groups(
  const cvm::real* group1_pos, const cvm::real* group2_pos,
  unsigned int numAtoms1, unsigned int numAtoms2,
  int en, int ed,
  const cvm::rvector inv_r0_vec,
  const cvm::rvector inv_r0sq_vec,
#if 0
  // TODO: Wait for https://github.com/Colvars/colvars/pull/919
  const cvm::system_boundary_conditions bc,
#endif
  cvm::real* group1_grad, cvm::real* group2_grad,
  cvm::real pairlist_tol,
  cvm::real pairlist_tol_l2_max,
  bool* d_pairlist,
  unsigned int* d_tbcount,
  cvm::real* d_coordnum_tmp,
  cvm::real* h_coordnum_out,
  const int flags,
  cudaStream_t stream,
  colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  if (numAtoms1 < numAtoms2) {
    return calc_value_coordnum_two_groups(
      group2_pos, group1_pos, numAtoms2, numAtoms1,
      en, ed, inv_r0_vec, inv_r0sq_vec,
#if 0
      // TODO: Wait for https://github.com/Colvars/colvars/pull/919
      const cvm::system_boundary_conditions bc,
#endif
      group2_grad, group1_grad,
      pairlist_tol, pairlist_tol_l2_max,
      d_pairlist, d_tbcount, d_coordnum_tmp,
      h_coordnum_out, flags, stream, cvmodule
    );
  }
  const cvm::real* pos1x = group1_pos;
  const cvm::real* pos1y = pos1x + numAtoms1;
  const cvm::real* pos1z = pos1y + numAtoms1;
  const cvm::real* pos2x = group2_pos;
  const cvm::real* pos2y = pos2x + numAtoms2;
  const cvm::real* pos2z = pos2y + numAtoms2;
  cvm::real* grad1x = group1_grad;
  cvm::real* grad1y = grad1x + numAtoms1;
  cvm::real* grad1z = grad1y + numAtoms1;
  cvm::real* grad2x = group2_grad;
  cvm::real* grad2y = grad2x + numAtoms2;
  cvm::real* grad2z = grad2y + numAtoms2;
  void* args[] = {
    &pos1x, &pos1y, &pos1z,
    &pos2x, &pos2y, &pos2z,
    &numAtoms1, &numAtoms2, &en, &ed,
    const_cast<cvm::rvector*>(&inv_r0_vec),
    const_cast<cvm::rvector*>(&inv_r0sq_vec),
#if 0
    // TODO: Wait for https://github.com/Colvars/colvars/pull/919
    const_cast<cvm::system_boundary_conditions*>(&bc),
#endif
    &grad1x, &grad1y, &grad1z,
    &grad2x, &grad2y, &grad2z,
    &pairlist_tol, &pairlist_tol_l2_max,
    &d_pairlist, &d_tbcount,
    &d_coordnum_tmp, &h_coordnum_out
  };
  constexpr int numGroup2BatchesPerBlock = 2;
  constexpr int group2WorkSize = default_block_size / numGroup2BatchesPerBlock;
  const unsigned int numBlocks = (numAtoms1 + default_block_size - 1) / default_block_size;
  void* kernel = nullptr;
#define CASE(N) case N: kernel = \
  (void*)computeCoordinationNumberTwoGroupsCUDAKernel1< \
    0, 0, default_block_size, group2WorkSize, numGroup2BatchesPerBlock, N>; break
  switch (flags) {
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_rebuild_pairlist +
         colvar::coordnum::ef_use_internal_pbc);

    CASE(colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_rebuild_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    default: {
      return cvmodule->error("Unimplemented flags in calc_value_coordnum_two_groups: " +
        cvm::to_str(flags) + "\n");
    }
  }
#undef CASE
  const unsigned int sharedMemBytes =
    (flags & colvar::coordnum::ef_use_pairlist) ?
    default_block_size * default_block_size * sizeof(bool) : 0;
  error_code |= checkGPUError(cudaLaunchKernel(
    kernel, dim3(numBlocks, 1, 1), dim3(default_block_size, 1, 1), args, sharedMemBytes, stream));
  return error_code;
}

template <int N, int M, int blockSize, int flags>
__global__ void computeCoordinationNumberGroupToCenterKernel(
  const cvm::real* __restrict posx,
  const cvm::real* __restrict posy,
  const cvm::real* __restrict posz,
  const cvm::rvector* __restrict com,
  const unsigned int numAtoms,
  const int en, const int ed,
  const cvm::rvector inv_r0_vec,
  const cvm::rvector inv_r0sq_vec,
#if 0
  // TODO: Wait for https://github.com/Colvars/colvars/pull/919
  const cvm::system_boundary_conditions bc,
#endif
  cvm::real* __restrict gx1,
  cvm::real* __restrict gy1,
  cvm::real* __restrict gz1,
  const cvm::real pairlist_tol,
  const cvm::real pairlist_tol_l2_max,
  bool* __restrict pairlist,
  unsigned int* __restrict tbcount,
  cvm::rvector* __restrict com_grad_tmp,
  cvm::rvector* __restrict com_grad_out,
  cvm::real* __restrict coordnum_tmp,
  cvm::real* __restrict coordnum_out) {
  constexpr const bool use_pairlist = flags & colvar::coordnum::ef_use_pairlist;
  constexpr const bool rebuild_pairlist = flags & colvar::coordnum::ef_rebuild_pairlist;
  constexpr const bool gradients = flags & colvar::coordnum::ef_gradients;
  constexpr const bool use_internal_pbc = flags & colvar::coordnum::ef_use_internal_pbc;
  static_assert(use_internal_pbc == true, "The CUDA kernel requires internal PBC.");
  __shared__ bool isLastBlockDone;
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  // Total coordnum
  cvm::real ei = 0;
  cvm::real com_grad_x, com_grad_y, com_grad_z;
  if constexpr (gradients) {
    com_grad_x = 0;
    com_grad_y = 0;
    com_grad_z = 0;
  }
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int gridSize = blockDim.x * gridDim.x;
  const cvm::real x2 = com->x;
  const cvm::real y2 = com->y;
  const cvm::real z2 = com->z;
  while (i < numAtoms) {
    cvm::real grad_x, grad_y, grad_z;
    if constexpr (gradients) {
      grad_x = 0;
      grad_y = 0;
      grad_z = 0;
    }
    const cvm::real x1 = posx[i];
    const cvm::real y1 = posy[i];
    const cvm::real z1 = posz[i];
#if 0
// TODO: Wait for https://github.com/Colvars/colvars/pull/919
    if constexpr (!use_pairlist) {
      ei += compute_pair_coordnum<flags>(
        inv_r0_vec, inv_r0sq_vec, en, ed,
        x1, y1, z1, x2, y2, z2,
        grad_x, grad_y, grad_z,
        com_grad_x,
        com_grad_y,
        com_grad_z,
        pairlist_tol, pairlist_tol_l2_max, bc);
    } else {
      if constexpr (!rebuild_pairlist) {
        const bool within = pairlist[i];
        if (within) {
          ei += compute_pair_coordnum<flags>(
            inv_r0_vec, inv_r0sq_vec, en, ed,
            x1, y1, z1, x2, y2, z2,
            grad_x, grad_y, grad_z,
            com_grad_x,
            com_grad_y,
            com_grad_z,
            pairlist_tol, pairlist_tol_l2_max, bc);
        }
      } else {
        const double f = compute_pair_coordnum<flags>(
          inv_r0_vec, inv_r0sq_vec, en, ed,
          x1, y1, z1, x2, y2, z2,
          grad_x, grad_y, grad_z,
          com_grad_x,
          com_grad_y,
          com_grad_z,
          pairlist_tol, pairlist_tol_l2_max, bc);
        pairlist[i] = f > 0.0;
        ei += f;
      }
    }
#endif
    if constexpr (gradients) {
      atomicAdd(&gx1[i], grad_x);
      atomicAdd(&gy1[i], grad_y);
      atomicAdd(&gz1[i], grad_z);
    }
    i += gridSize;
  }
  typedef cub::BlockReduce<cvm::real, blockSize> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  const cvm::real total_e = BlockReduce(temp_storage).Sum(ei); __syncthreads();
  const cvm::real total_com_grad_x = BlockReduce(temp_storage).Sum(com_grad_x); __syncthreads();
  const cvm::real total_com_grad_y = BlockReduce(temp_storage).Sum(com_grad_y); __syncthreads();
  const cvm::real total_com_grad_z = BlockReduce(temp_storage).Sum(com_grad_z); __syncthreads();
  if (threadIdx.x == 0) {
    atomicAdd(coordnum_tmp, total_e);
    atomicAdd(&com_grad_tmp->x, total_com_grad_x);
    atomicAdd(&com_grad_tmp->y, total_com_grad_y);
    atomicAdd(&com_grad_tmp->z, total_com_grad_z);
    __threadfence();
    unsigned int value = atomicInc(tbcount, gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      *coordnum_out = *coordnum_tmp;
      *com_grad_out = *com_grad_tmp;
      *coordnum_tmp = 0;
      com_grad_tmp->x = 0;
      com_grad_tmp->y = 0;
      com_grad_tmp->z = 0;
      tbcount[0] = 0;
    }
  }
}

int calc_value_coordnum_group_to_com(
  const cvm::real* group_pos,
  const cvm::rvector* com,
  unsigned int numAtoms,
  int en, int ed,
  const cvm::rvector inv_r0_vec,
  const cvm::rvector inv_r0sq_vec,
#if 0
  // TODO: Wait for https://github.com/Colvars/colvars/pull/919
  const cvm::system_boundary_conditions bc,
#endif
  cvm::real* group_grad,
  cvm::real pairlist_tol,
  cvm::real pairlist_tol_l2_max,
  bool* d_pairlist,
  unsigned int* d_tbcount,
  cvm::rvector* d_com_grad_tmp,
  cvm::rvector* d_com_grad_out,
  cvm::real* d_coordnum_tmp,
  cvm::real* h_coordnum_out,
  const int flags,
  cudaStream_t stream,
  colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  const cvm::real* posx = group_pos;
  const cvm::real* posy = posx+ numAtoms;
  const cvm::real* posz = posy + numAtoms;
  cvm::real* gradx = group_grad;
  cvm::real* grady = gradx + numAtoms;
  cvm::real* gradz = grady + numAtoms;
  void* args[] = {
    &posx, &posy, &posz,
    &com, &numAtoms, &en, &ed,
    const_cast<cvm::rvector*>(&inv_r0_vec),
    const_cast<cvm::rvector*>(&inv_r0sq_vec),
#if 0
    // TODO: Wait for https://github.com/Colvars/colvars/pull/919
    const_cast<cvm::system_boundary_conditions*>(&bc),
#endif
    &gradx, &grady, &gradz,
    &pairlist_tol, &pairlist_tol_l2_max,
    &d_pairlist, &d_tbcount,
    &d_com_grad_tmp, &d_com_grad_out,
    &d_coordnum_tmp, &h_coordnum_out};
  const unsigned int numBlocks = (numAtoms + default_block_size - 1) / default_block_size;
  void* kernel = nullptr;
#define CASE(N) case N: kernel = \
  (void*)computeCoordinationNumberGroupToCenterKernel< \
    0, 0, default_block_size, N>; break
  switch (flags) {
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_rebuild_pairlist +
         colvar::coordnum::ef_use_internal_pbc);

    CASE(colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_rebuild_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    default: {
      return cvmodule->error("Unimplemented flags in calc_value_coordnum_group_to_com: " +
        cvm::to_str(flags) + "\n");
    }
  }
#undef CASE
  const unsigned int sharedMemBytes = 0;
  error_code |= checkGPUError(cudaLaunchKernel(
    kernel, dim3(numBlocks, 1, 1), dim3(default_block_size, 1, 1), args, sharedMemBytes, stream));
  return error_code;
}

template <int N, int M, int flags>
void computeCoordinationNumberGroupTwoCOMsKernel(
  const cvm::rvector* __restrict com1,
  const cvm::rvector* __restrict com2,
  int en, int ed,
  const cvm::rvector inv_r0_vec,
  const cvm::rvector inv_r0sq_vec,
#if 0
  // TODO: Wait for https://github.com/Colvars/colvars/pull/919
  const cvm::system_boundary_conditions bc,
#endif
  cvm::real pairlist_tol,
  cvm::real pairlist_tol_l2_max,
  bool* __restrict pairlist,
  cvm::rvector* __restrict com1_grad_out,
  cvm::rvector* __restrict com2_grad_out,
  cvm::real* __restrict h_coordnum_out) {
  constexpr const bool use_pairlist = flags & colvar::coordnum::ef_use_pairlist;
  constexpr const bool rebuild_pairlist = flags & colvar::coordnum::ef_rebuild_pairlist;
  constexpr const bool gradients = flags & colvar::coordnum::ef_gradients;
  constexpr const bool use_internal_pbc = flags & colvar::coordnum::ef_use_internal_pbc;
  static_assert(use_internal_pbc == true, "The CUDA kernel requires internal PBC.");
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  // Use only 1 thread
  if (i < 1) {
    const cvm::real x1 = com1->x;
    const cvm::real y1 = com1->y;
    const cvm::real z1 = com1->z;
    const cvm::real x2 = com2->x;
    const cvm::real y2 = com2->y;
    const cvm::real z2 = com2->z;
    cvm::real com1_grad_x;
    cvm::real com1_grad_y;
    cvm::real com1_grad_z;
    cvm::real com2_grad_x;
    cvm::real com2_grad_y;
    cvm::real com2_grad_z;
    if constexpr (gradients) {
      com1_grad_x = 0;
      com1_grad_y = 0;
      com1_grad_z = 0;
      com2_grad_x = 0;
      com2_grad_y = 0;
      com2_grad_z = 0;
    }
#if 0
// TODO: Wait for https://github.com/Colvars/colvars/pull/919
    if constexpr (!use_pairlist) {
      (*h_coordnum_out) = compute_pair_coordnum<flags>(
        inv_r0_vec, inv_r0sq_vec, en, ed,
        x1, y1, z1, x2, y2, z2,
        com1_grad_x, com1_grad_y, com1_grad_z,
        com2_grad_x, com2_grad_y, com2_grad_z,
        pairlist_tol, pairlist_tol_l2_max, bc);
    } else {
      if constexpr (!rebuild_pairlist) {
        const bool within = pairlist[i];
        if (within) {
          (*h_coordnum_out) = compute_pair_coordnum<flags>(
            inv_r0_vec, inv_r0sq_vec, en, ed,
            x1, y1, z1, x2, y2, z2,
            com1_grad_x, com1_grad_y, com1_grad_z,
            com2_grad_x, com2_grad_y, com2_grad_z,
            pairlist_tol, pairlist_tol_l2_max, bc);
        }
      } else {
        const double f = compute_pair_coordnum<flags>(
            inv_r0_vec, inv_r0sq_vec, en, ed,
            x1, y1, z1, x2, y2, z2,
            com1_grad_x, com1_grad_y, com1_grad_z,
            com2_grad_x, com2_grad_y, com2_grad_z,
            pairlist_tol, pairlist_tol_l2_max, bc);
        pairlist[i] = f > 0.0;
        (*h_coordnum_out) = f;
      }
    }
#endif
  }
}

int calc_value_coordnum_com_to_com(
  const cvm::rvector* d_com1,
  const cvm::rvector* d_com2,
  int en, int ed,
  const cvm::rvector inv_r0_vec,
  const cvm::rvector inv_r0sq_vec,
#if 0
  // TODO: Wait for https://github.com/Colvars/colvars/pull/919
  const cvm::system_boundary_conditions bc,
#endif
  cvm::real pairlist_tol,
  cvm::real pairlist_tol_l2_max,
  bool* d_pairlist,
  cvm::rvector* d_com1_grad_out,
  cvm::rvector* d_com2_grad_out,
  cvm::real* h_coordnum_out,
  const int flags,
  cudaStream_t stream,
  colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  void* args[] = {
    &d_com1, &d_com2, &en, &ed,
    const_cast<cvm::rvector*>(&inv_r0_vec),
    const_cast<cvm::rvector*>(&inv_r0sq_vec),
#if 0
    // TODO: Wait for https://github.com/Colvars/colvars/pull/919
    const_cast<cvm::system_boundary_conditions*>(&bc),
#endif
    &pairlist_tol, &pairlist_tol_l2_max,
    &d_pairlist,
    &d_com1_grad_out,
    &d_com2_grad_out,
    &h_coordnum_out};
  void* kernel = nullptr;
#define CASE(N) case N: kernel = \
  (void*)computeCoordinationNumberGroupTwoCOMsKernel< \
    0, 0, N>; break
  switch (flags) {
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_rebuild_pairlist +
         colvar::coordnum::ef_use_internal_pbc);

    CASE(colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_rebuild_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    default: {
      return cvmodule->error("Unimplemented flags in calc_value_coordnum_com_to_com: " +
        cvm::to_str(flags) + "\n");
    }
  }
#undef CASE
  const unsigned int sharedMemBytes = 0;
  error_code |= checkGPUError(cudaLaunchKernel(
    kernel, dim3(1, 1, 1), dim3(default_block_size, 1, 1), args, sharedMemBytes, stream));
  return error_code;
}

__inline__ __device__ size_t computeGlobalPairlistIDSelfGroup(
  size_t iid_global, size_t jid_global, size_t numAtoms) {
  const size_t iid = min(iid_global, jid_global);
  const size_t jid = max(iid_global, jid_global);
  return iid * (2 * numAtoms - 1 - iid) / 2 + (jid - iid - 1);
}

// WARNING: To match the CPU data layout, the GPU pairlist implementation is very inefficient!
template <int N, int M, int blockSize, int tileSize, int flags>
__global__ void computeCoordinationNumberSelfGroupCUDAKernel1(
  const cvm::real* __restrict pos1x,
  const cvm::real* __restrict pos1y,
  const cvm::real* __restrict pos1z,
  const unsigned int numAtoms1,
  const int en, const int ed,
  const cvm::rvector inv_r0_vec,
  const cvm::rvector inv_r0sq_vec,
#if 0
  // TODO: Wait for https://github.com/Colvars/colvars/pull/919
  const cvm::system_boundary_conditions bc,
#endif
  cvm::real* __restrict gx1,
  cvm::real* __restrict gy1,
  cvm::real* __restrict gz1,
  const unsigned int* __restrict tilesList,
  const unsigned int* __restrict tilesListStart,
  const unsigned int* __restrict tilesListSizes,
  const cvm::real pairlist_tol,
  const cvm::real pairlist_tol_l2_max,
  bool* __restrict pairlist,
  unsigned int* __restrict tbcount,
  cvm::real* __restrict coordnum_tmp,
  cvm::real* __restrict coordnum_out) {
#if 0
  // TODO: Optimize for static en and ed. Wait for https://github.com/Colvars/colvars/pull/926.
  constexpr const bool static_exponents = (N > 0) && (M > 0);
#endif
  constexpr const bool use_pairlist = flags & colvar::coordnum::ef_use_pairlist;
  constexpr const bool rebuild_pairlist = flags & colvar::coordnum::ef_rebuild_pairlist;
  constexpr const bool gradients = flags & colvar::coordnum::ef_gradients;
  constexpr const bool use_internal_pbc = flags & colvar::coordnum::ef_use_internal_pbc;
  constexpr unsigned int numTilesPerBlock = blockSize / tileSize;
  __shared__ double3 shJGrad[numTilesPerBlock][tileSize];
  extern __shared__ unsigned int globalJIDs_buffer[];
  unsigned int (&globalJIDs)[numTilesPerBlock][tileSize] =
    *reinterpret_cast<unsigned int (*)[numTilesPerBlock][tileSize]>(globalJIDs_buffer);
  __shared__ bool isLastBlockDone;
  cvm::real coordnum_tb = 0;
  static constexpr const unsigned int half_tile_size = tileSize / 2;
  // Number of blocks required to iterate over group1
  const unsigned int numBlocksInGroup1 = (numAtoms1 + blockSize - 1) / blockSize;
  const unsigned int numTilesInGroup1 = (numAtoms1 + tileSize - 1) / tileSize;
  namespace cg = cooperative_groups;
  const cg::thread_block thb = cg::this_thread_block();
  const auto tilePartition = cg::tiled_partition<tileSize>(thb);
  // Which tile is the current thread working on?
  const unsigned int tileIndexInBlock = tilePartition.meta_group_rank();
  // The thread id in a tile
  const unsigned int threadIndexInTile = tilePartition.thread_rank();
  for (unsigned int i = blockIdx.x; i < numBlocksInGroup1; i += gridDim.x) {
    unsigned int iTileIndexInGrid = i * numTilesPerBlock + tileIndexInBlock;
    const bool isItileValid = iTileIndexInGrid < numTilesInGroup1;
    if (isItileValid) {
      const unsigned int tid = i * blockDim.x + threadIdx.x;
      const bool mask_i = tid < numAtoms1;
      const cvm::real x1 = mask_i ? pos1x[tid] : 0;
      const cvm::real y1 = mask_i ? pos1y[tid] : 0;
      const cvm::real z1 = mask_i ? pos1z[tid] : 0;
      double3 iGrad{0, 0, 0};
      if constexpr (use_pairlist) {
        globalJIDs[tileIndexInBlock][threadIndexInTile] = tid;
      }
      // Self tile
      shJGrad[tileIndexInBlock][threadIndexInTile].x = 0;
      shJGrad[tileIndexInBlock][threadIndexInTile].y = 0;
      shJGrad[tileIndexInBlock][threadIndexInTile].z = 0;
      tilePartition.sync();
      // Self tiles
      #pragma unroll
      for (unsigned int t = 1; t < half_tile_size; ++t) {
        // NAMD/OpenMM style swizzling
        const unsigned int jid = (t + threadIndexInTile) & (tileSize - 1);
        const bool mask_t = tilePartition.shfl(mask_i, jid);
        unsigned int pairlistID;
        bool pairlist_elem;
        unsigned int jid_global;
        const cvm::real x2 = tilePartition.shfl(x1, jid);
        const cvm::real y2 = tilePartition.shfl(y1, jid);
        const cvm::real z2 = tilePartition.shfl(z1, jid);
        if (mask_i && mask_t) {
          if constexpr (use_pairlist) {
            jid_global = globalJIDs[tileIndexInBlock][jid];
            pairlistID = computeGlobalPairlistIDSelfGroup(tid, jid_global, numAtoms1);
          }
          if constexpr (use_pairlist && !rebuild_pairlist) {
            pairlist_elem = pairlist[pairlistID];
          }
#if 0
          // TODO: Wait for https://github.com/Colvars/colvars/pull/919
          const auto partial = compute_pair_coordnum<flags>(
            inv_r0_vec, inv_r0sq_vec, en, ed,
            x1, y1, z1, x2, y2, z2,
            iGrad.x, iGrad.y, iGrad.z,
            shJGrad[tileIndexInBlock][jid].x,
            shJGrad[tileIndexInBlock][jid].y,
            shJGrad[tileIndexInBlock][jid].z,
            pairlist_tol, pairlist_tol_l2_max, bc);
          coordnum_tb += partial;
          if constexpr (use_pairlist && rebuild_pairlist) {
            pairlist_elem = partial > 0.0 ? true : false;
          }
#endif
          if constexpr (use_pairlist && rebuild_pairlist) {
            pairlist[pairlistID] = pairlist_elem;
          }
        }
        tilePartition.sync();
      }
      // Last loop: t == block_size / 2
      {
        const unsigned int jid = (half_tile_size + threadIndexInTile) & (tileSize - 1);
        unsigned int pairlistID;
        bool pairlist_elem;
        const bool mask_t = tilePartition.shfl(mask_i, jid);
        const cvm::real x2 = tilePartition.shfl(x1, jid);
        const cvm::real y2 = tilePartition.shfl(y1, jid);
        const cvm::real z2 = tilePartition.shfl(z1, jid);
        if (jid > threadIndexInTile) {
          if (mask_i && mask_t) {
            if constexpr (use_pairlist) {
              const unsigned int jid_global = globalJIDs[tileIndexInBlock][jid];
              pairlistID = computeGlobalPairlistIDSelfGroup(tid, jid_global, numAtoms1);
            }
            if constexpr (use_pairlist && !rebuild_pairlist) {
              pairlist_elem = pairlist[pairlistID];
            }
#if 0
            // TODO: Wait for https://github.com/Colvars/colvars/pull/919
            const auto partial = compute_pair_coordnum<flags>(
              inv_r0_vec, inv_r0sq_vec, en, ed,
              x1, y1, z1, x2, y2, z2,
              iGrad.x, iGrad.y, iGrad.z,
              shJGrad[tileIndexInBlock][jid].x,
              shJGrad[tileIndexInBlock][jid].y,
              shJGrad[tileIndexInBlock][jid].z,
              pairlist_tol, pairlist_tol_l2_max, bc);
            coordnum_tb += partial;
            if constexpr (use_pairlist && rebuild_pairlist) {
              pairlist_elem = partial > 0.0 ? true : false;
            }
#endif
            if constexpr (use_pairlist && rebuild_pairlist) {
              pairlist[pairlistID] = pairlist_elem;
            }
          }
        }
        tilePartition.sync();
      }
      if (mask_i) {
        atomicAdd(&gx1[tid], shJGrad[tileIndexInBlock][threadIndexInTile].x);
        atomicAdd(&gy1[tid], shJGrad[tileIndexInBlock][threadIndexInTile].y);
        atomicAdd(&gz1[tid], shJGrad[tileIndexInBlock][threadIndexInTile].z);
      }
      // Iterate over other tiles
      const unsigned int jTileStart = tilesListStart[iTileIndexInGrid];
      const unsigned int numJTiles = tilesListSizes[iTileIndexInGrid];
      const unsigned int jTileEnd = jTileStart + numJTiles;
      for (unsigned int l = jTileStart; l < jTileEnd; ++l) {
        const unsigned int jTileIndex = tilesList[l];
        // Fetch atom j from i-tile
        const unsigned int jid_global = jTileIndex * tileSize + threadIndexInTile;
        const bool mask_j = jid_global < numAtoms1;
        const cvm::real jx2 = mask_j ? pos1x[jid_global] : 0;
        const cvm::real jy2 = mask_j ? pos1y[jid_global] : 0;
        const cvm::real jz2 = mask_j ? pos1z[jid_global] : 0;
        // Reset the gradients
        shJGrad[tileIndexInBlock][threadIndexInTile].x = 0;
        shJGrad[tileIndexInBlock][threadIndexInTile].y = 0;
        shJGrad[tileIndexInBlock][threadIndexInTile].z = 0;
        if constexpr (use_pairlist) {
          globalJIDs[tileIndexInBlock][threadIndexInTile] = jid_global;
        }
        tilePartition.sync();
        #pragma unroll
        for (unsigned int t = 0; t < tileSize; ++t) {
          const unsigned int jid = t ^ threadIndexInTile;
          const bool mask_t = tilePartition.shfl(mask_j, jid);
          unsigned int pairlistID;
          bool pairlist_elem;
          const cvm::real x2 = tilePartition.shfl(jx2, jid);
          const cvm::real y2 = tilePartition.shfl(jy2, jid);
          const cvm::real z2 = tilePartition.shfl(jz2, jid);
          if (mask_i && mask_t) {
            if constexpr (use_pairlist) {
              pairlistID = computeGlobalPairlistIDSelfGroup(
                tid, globalJIDs[tileIndexInBlock][jid], numAtoms1);
            }
            if constexpr (use_pairlist && !rebuild_pairlist) {
              pairlist_elem = pairlist[pairlistID];
            }
#if 0
            // TODO: Wait for https://github.com/Colvars/colvars/pull/919
            const auto partial = compute_pair_coordnum<flags>(
              inv_r0_vec, inv_r0sq_vec, en, ed,
              x1, y1, z1, x2, y2, z2,
              iGrad.x, iGrad.y, iGrad.z,
              shJGrad[tileIndexInBlock][jid].x,
              shJGrad[tileIndexInBlock][jid].y,
              shJGrad[tileIndexInBlock][jid].z,
              pairlist_tol, pairlist_tol_l2_max, bc);
            coordnum_tb += partial;
            if constexpr (use_pairlist && rebuild_pairlist) {
              pairlist_elem = partial > 0.0 ? true : false;
            }
#endif
            if constexpr (use_pairlist && rebuild_pairlist) {
              pairlist[pairlistID] = pairlist_elem;
            }
          }
          tilePartition.sync();
        }
        if (mask_j) {
          atomicAdd(&gx1[jid_global], shJGrad[tileIndexInBlock][threadIndexInTile].x);
          atomicAdd(&gy1[jid_global], shJGrad[tileIndexInBlock][threadIndexInTile].y);
          atomicAdd(&gz1[jid_global], shJGrad[tileIndexInBlock][threadIndexInTile].z);
        }
      }
      if (mask_i) {
        atomicAdd(&gx1[tid], iGrad.x);
        atomicAdd(&gy1[tid], iGrad.y);
        atomicAdd(&gz1[tid], iGrad.z);
      }
    }
  }
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  // Reduction for energy
  typedef cub::BlockReduce<cvm::real, blockSize> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  coordnum_tb = BlockReduce(temp_storage).Sum(coordnum_tb); __syncthreads();
  if (threadIdx.x == 0) {
    atomicAdd(coordnum_tmp, coordnum_tb);
    __threadfence();
    unsigned int value = atomicInc(tbcount, gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      *coordnum_out = *coordnum_tmp;
      *coordnum_tmp = 0;
      tbcount[0] = 0;
    }
  }
}

int calc_value_coordnum_self_group(
  const cvm::real* group_pos,
  unsigned int numAtoms,
  int en, int ed,
  const cvm::rvector inv_r0_vec,
  const cvm::rvector inv_r0sq_vec,
#if 0
  // TODO: Wait for https://github.com/Colvars/colvars/pull/919
  const cvm::system_boundary_conditions bc,
#endif
  cvm::real* group_grad,
  const unsigned int* d_tilesList,
  const unsigned int* d_tilesListStart,
  const unsigned int* d_tilesListSizes,
  cvm::real pairlist_tol,
  cvm::real pairlist_tol_l2_max,
  bool* d_pairlist,
  unsigned int* d_tbcount,
  cvm::real* d_coordnum_tmp,
  cvm::real* h_coordnum_out,
  const int flags,
  cudaStream_t stream,
  colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  const cvm::real* posx = group_pos;
  const cvm::real* posy = posx+ numAtoms;
  const cvm::real* posz = posy + numAtoms;
  cvm::real* gradx = group_grad;
  cvm::real* grady = gradx + numAtoms;
  cvm::real* gradz = grady + numAtoms;
  void* args[] = {
    &posx, &posy, &posz,
    &numAtoms, &en, &ed,
    const_cast<cvm::rvector*>(&inv_r0_vec),
    const_cast<cvm::rvector*>(&inv_r0sq_vec),
#if 0
    // TODO: Wait for https://github.com/Colvars/colvars/pull/919
    const_cast<cvm::system_boundary_conditions*>(&bc),
#endif
    &gradx, &grady, &gradz,
    const_cast<unsigned int**>(&d_tilesList),
    const_cast<unsigned int**>(&d_tilesListStart),
    const_cast<unsigned int**>(&d_tilesListSizes),
    &pairlist_tol, &pairlist_tol_l2_max,
    &d_pairlist, &d_tbcount,
    &d_coordnum_tmp, &h_coordnum_out
  };
  const unsigned int numBlocks = (numAtoms + default_block_size - 1) / default_block_size;
  void* kernel = nullptr;
  const int gpu_warp_size = cvmodule->proxy->gpu_warp_size();
  // NOTE: For CUDA, we only support 32 as a warp size, but for HIP,
  // there could be two different warp sizes (32 and 64).
#if defined (COLVARS_CUDA)
#define CASE(N) case N: { \
  switch (gpu_warp_size) { \
    case 32: kernel = (void*)computeCoordinationNumberSelfGroupCUDAKernel1<0, 0, default_block_size, 32, N>; break;\
    default: return cvmodule->error("Unsupported warp size in calc_value_coordnum_self_group: " + cvm::to_str(gpu_warp_size) + "\n", COLVARS_BUG_ERROR);\
  } \
  } break
#elif defined (COLVARS_HIP)
#define CASE(N) case N: { \
  switch (gpu_warp_size) { \
    case 32: kernel = (void*)computeCoordinationNumberSelfGroupCUDAKernel1<0, 0, default_block_size, 32, N>; break;\
    case 64: kernel = (void*)computeCoordinationNumberSelfGroupCUDAKernel1<0, 0, default_block_size, 64, N>; break;\
    default: return cvmodule->error("Unsupported warp size in calc_value_coordnum_self_group: " + cvm::to_str(gpu_warp_size) + "\n", COLVARS_BUG_ERROR);\
  } \
  } break
#endif
  switch (flags) {
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_gradients +
         colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_rebuild_pairlist +
         colvar::coordnum::ef_use_internal_pbc);

    CASE(colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    CASE(colvar::coordnum::ef_use_pairlist +
         colvar::coordnum::ef_rebuild_pairlist +
         colvar::coordnum::ef_use_internal_pbc);
    default: {
      return cvmodule->error("Unimplemented flags in calc_value_coordnum_self_group: " +
        cvm::to_str(flags) + "\n");
    }
  }
#undef CASE
  const unsigned int sharedMemBytes = (flags & colvar::coordnum::ef_use_pairlist) ? default_block_size * sizeof(unsigned int) : 0;
  error_code |= checkGPUError(cudaLaunchKernel(
    kernel, dim3(numBlocks, 1, 1), dim3(default_block_size, 1, 1), args, sharedMemBytes, stream));
  return error_code;
}

}
#endif // defined(COLVARS_CUDA) || defined (COLVARS_HIP)
