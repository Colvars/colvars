#include "colvar_rotation_derivative.h"
#include "colvarproxy.h"

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

#include "cuda/colvar_rotation_derivative_kernel.h"

namespace colvars_gpu {
rotation_derivative_gpu::rotation_derivative_gpu():
  m_rot(nullptr), m_d_pos1(nullptr), m_d_pos2(nullptr),
  m_num_atoms_pos1(0),
  m_num_atoms_pos2(0),
  tmp_Q0Q0(nullptr), tmp_Q0Q0_L(nullptr) {
}

int rotation_derivative_gpu::init(
  const colvars_gpu::rotation_gpu* rot,
  const cvm::real* d_pos1,
  const cvm::real* d_pos2,
  const size_t num_atoms_pos1,
  const size_t num_atoms_pos2) {
  int error_code = COLVARS_OK;
  m_rot = rot;
  m_d_pos1 = d_pos1;
  m_d_pos2 = d_pos2;
  m_num_atoms_pos1 = num_atoms_pos1;
  m_num_atoms_pos2 = num_atoms_pos2;
  colvarproxy* p = cvm::main()->proxy;
  error_code |= p->reallocate_device(&tmp_Q0Q0, 4 * 4);
  error_code |= p->reallocate_device(&tmp_Q0Q0_L, 4 * 4 * 4);
  return error_code;
}

rotation_derivative_gpu::~rotation_derivative_gpu() {
  colvarproxy* p = cvm::main()->proxy;
  p->deallocate_device(&tmp_Q0Q0);
  p->deallocate_device(&tmp_Q0Q0_L);
  tmp_Q0Q0 = nullptr;
  tmp_Q0Q0_L = nullptr;
}

int rotation_derivative_gpu::add_prepare_derivative_nodes(
  rotation_derivative_dldq require_dl_dq,
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  cudaGraphNode_t node;
  error_code |= colvars_gpu::prepare_derivative(
    require_dl_dq, m_rot->d_S_eigval,
    m_rot->d_S_eigvec, tmp_Q0Q0,
    tmp_Q0Q0_L, node, graph, {});
  nodes_map["prepare_rotation_derivative"] = node;
  return error_code;
}

}

#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
