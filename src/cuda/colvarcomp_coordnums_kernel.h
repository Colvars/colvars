#ifndef COLVARCOMP_COORDNUMS_KERNEL_H
#define COLVARCOMP_COORDNUMS_KERNEL_H

#include "colvarmodule.h"
#include "colvar_gpu_support.h"

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

namespace colvars_gpu {
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
  colvarmodule* cvmodule);

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
  colvarmodule* cvmodule);

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
  colvarmodule* cvmodule);

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
  colvarmodule* cvmodule);
}

#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
#endif // COLVARCOMP_COORDNUMS_KERNEL_H
