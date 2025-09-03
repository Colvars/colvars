#ifndef COLVARATOMS_GPU_H
#define COLVARATOMS_GPU_H

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvar_rotation_derivative.h"

namespace colvars_gpu {

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
struct colvaratoms_gpu_buffer_t {
  /// \brief GPU atom proxy indices (size: num_atoms)
  int* d_atoms_index;
  /// \brief GPU atom positions (size: 3 * num_atoms)
  cvm::real* d_atoms_pos;
  /// \brief GPU atom charges (size: num_atoms)
  cvm::real* d_atoms_charge;
  /// \brief GPU atom velocities (size: 3 * num_atoms)
  cvm::real* d_atoms_vel;
  /// \brief GPU atom mass (size: num_atoms)
  cvm::real* d_atoms_mass;
  /// \brief GPU atom gradients (size: 3 * num_atoms)
  cvm::real* d_atoms_grad;
  /// \brief GPU atom total forces (size: 3 * num_atoms)
  cvm::real* d_atoms_total_force;
  /// \brief Atom masses divided by total mass (size: num_atoms)
  cvm::real* d_atoms_weight;
  /// \brief GPU atom applied force
  cvm::real* d_atoms_applied_force;
  /// \brief GPU fit gradients
  cvm::real* d_fit_gradients;
  /// \brief GPU reference coordinates for f_ag_center or f_ag_rotate
  cvm::real* d_ref_pos;
  /// \brief GPU atom positions (size: 3 * num_atoms)
  cvm::real* d_atoms_pos_unrotated;
  /// \brief GPU center-of-mass
  cvm::rvector* d_com;
  /// \brief GPU temporary buffer for COM, used for avoiding memset
  cvm::rvector* d_com_tmp;
  /// \brief GPU center-of-geometry
  cvm::rvector* d_cog;
  /// \brief GPU temporary buffer for COG, used for avoiding memset
  cvm::rvector* d_cog_tmp;
  /// \brief GPU center of geometry before any fitting
  cvm::rvector* d_cog_orig;
  /// \brief GPU atomic counter for reduction
  unsigned int* d_com_cog_tbcount;
  /// \brief Center-of-mass on the host-pinned memory for CPU compatibility
  cvm::rvector* h_com;
  /// \brief Center-of-geometry on the host-pinned memory for CPU compatibility
  cvm::rvector* h_cog;
  /// \brief Center-of-geometry before any fitting on the host-pinned memory for CPU compatibility
  cvm::rvector* h_cog_orig;
  /// \brief GPU center of geometry of the reference coordinates
  cvm::rvector* d_ref_pos_cog;
};

struct colvaratoms_gpu_calc_fit_info_t {
  double3* d_atom_grad;
  double4* d_sum_dxdq;
  unsigned int* d_tbcount;
};

struct colvaratoms_gpu_debug_graph_t {
  bool initialized;
  cudaGraph_t graph_calc_required_properties;
  cudaGraphExec_t graph_exec_calc_required_properties;
};

class colvaratoms_gpu {
public:
  colvaratoms_gpu();
  ~colvaratoms_gpu();
  int init_gpu();
  int destroy_gpu();
  int sync_to_gpu_buffers(const cvm::atom_group* cpu_atoms);
  int clear_gpu_buffers(const cvm::atom_group* cpu_atoms);
  int add_reset_atoms_data_nodes(
    const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map);
  int add_read_positions_nodes(
    const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map);
  int add_calc_required_properties_nodes(
    const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
    const std::vector<cudaGraphNode_t>& extra_initial_dependencies = {});
  int add_update_cpu_buffers_nodes(
    cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map);
  int after_read_data_sync(
    cvm::atom_group* cpu_atoms, bool copy_to_cpu, cudaStream_t stream);
  int begin_apply_force_gpu();
  int add_apply_force_nodes(
    const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
    const std::vector<cudaGraphNode_t>& extra_initial_dependencies = {});
  int add_calc_fit_gradients_nodes(
    cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
    bool use_cpu_buffers = false);
  // For debug gradients
  int read_positions_gpu_debug(
    cvm::atom_group* cpu_atoms, bool change_fitting_group, size_t change_atom_i,
    int xyz, bool to_cpu, double sign, cudaStream_t stream);
  int calc_required_properties_gpu_debug(
    cvm::atom_group* cpu_atoms, bool to_cpu, cudaStream_t stream);
  void do_feature_side_effects_gpu(
    cvm::atom_group* cpu_atoms, int id);
  int setup_rotation(const cvm::atom_group* cpu_atoms);
  int setup_rotation_derivative(const cvm::atom_group* cpu_atoms);
  colvaratoms_gpu_buffer_t& get_gpu_buffers() { return gpu_buffers; }
  int read_total_forces(cvm::atom_group* cpu_atoms);
  void apply_colvar_force_from_cpu(cvm::real const& cpu_force);
  void set_use_cpu_group_force(bool yesno) { use_group_force = yesno; }
private:
  colvaratoms_gpu_buffer_t gpu_buffers;
  /// \brief Temporary variables for calc_fit_gradients GPU kernel
  colvaratoms_gpu_calc_fit_info_t calc_fit_gradients_gpu_info;
  colvaratoms_gpu_calc_fit_info_t calc_fit_forces_gpu_info;
  /// \brief Separate CUDA graphs for supporting debug gradients
  colvaratoms_gpu_debug_graph_t debug_graphs;
  /// \brief For intercepting the forces applied from the CPU interface
  cvm::real* h_sum_applied_colvar_force;
  /// \brief If the CPU code path use apply_colvar_force(),
  /// this will be set to true, and then reset to false in begin_apply_force_gpu()
  bool use_apply_colvar_force;
  /// \brief If the CPU code path use group_force_object,
  /// this will be set to true, and then reset to false in begin_apply_force_gpu()
  bool use_group_force;
  /// \brief GPU rotation object
  colvars_gpu::rotation_gpu rot_gpu;
  /// \brief GPU Rotation derivative;
  colvars_gpu::rotation_derivative_gpu* rot_deriv_gpu;
};

#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
}

#endif // COLVARATOMS_GPU_H
