#include "colvaratoms.h"
#include "colvar_gpu_support.h"
#include "colvardeps.h"
#include "colvarproxy.h"
#include "colvarmodule.h"
#include "colvarcomp.h"

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
#include "cuda/colvaratoms_kernel.h"
#endif

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
int cvm::atom_group::init_gpu() {
  int error_code = COLVARS_OK;
  colvarproxy *p = cvm::main()->proxy;
  // error_code |= checkGPUError(cudaStreamCreate(&stream));
  error_code |= p->reallocate_device(&gpu_buffers.d_com, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_cog, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_cog_orig, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_ref_pos_cog, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_com_cog_tbcount, 1);
  error_code |= p->clear_device_array(gpu_buffers.d_com_cog_tbcount, 1);
  error_code |= p->reallocate_host(&gpu_buffers.h_com, 1);
  error_code |= p->reallocate_host(&gpu_buffers.h_cog, 1);
  error_code |= p->reallocate_host(&gpu_buffers.h_cog_orig, 1);
  rot_deriv_gpu = nullptr;
  // error_code |= checkGPUError(cudaStreamCreate(&stream_ag_force));
  // std::memset(&calc_fit_gradients_info, 0, sizeof(calc_fit_gradients_info));
  error_code |= p->reallocate_device(&calc_fit_gradients_gpu_info.d_atom_grad, 1);
  error_code |= p->reallocate_device(&calc_fit_gradients_gpu_info.d_sum_dxdq, 1);
  error_code |= p->reallocate_device(&calc_fit_gradients_gpu_info.d_tbcount, 1);
  error_code |= p->clear_device_array(calc_fit_gradients_gpu_info.d_tbcount, 1);
  error_code |= p->reallocate_device(&calc_fit_forces_gpu_info.d_atom_grad, 1);
  error_code |= p->reallocate_device(&calc_fit_forces_gpu_info.d_sum_dxdq, 1);
  error_code |= p->reallocate_device(&calc_fit_forces_gpu_info.d_tbcount, 1);
  error_code |= p->clear_device_array(calc_fit_forces_gpu_info.d_tbcount, 1);
  error_code |= p->reallocate_host(&h_sum_applied_colvar_force, 1);
  use_apply_colvar_force = false;
  use_group_force = false;
  if (debug_graphs.graph_calc_required_properties) {
    error_code |= checkGPUError(cudaGraphDestroy(
      debug_graphs.graph_calc_required_properties));
  }
  if (debug_graphs.graph_exec_calc_required_properties) {
    error_code |= checkGPUError(cudaGraphExecDestroy(
      debug_graphs.graph_exec_calc_required_properties));
  }
  debug_graphs.initialized = false;
  return error_code;
}

int cvm::atom_group::destroy_gpu() {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvm::main()->proxy;
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_index);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_pos);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_charge);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_vel);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_mass);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_grad);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_total_force);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_weight);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_applied_force);
  error_code |= p->deallocate_device(&gpu_buffers.d_fit_gradients);
  error_code |= p->deallocate_device(&gpu_buffers.d_ref_pos);
  if (rot_deriv_gpu) {
    rot_deriv_gpu->~rotation_derivative_gpu();
    error_code |= p->deallocate_host(&rot_deriv_gpu);
  }
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_pos_unrotated);
  error_code |= p->deallocate_device(&gpu_buffers.d_com);
  error_code |= p->deallocate_device(&gpu_buffers.d_cog);
  error_code |= p->deallocate_device(&gpu_buffers.d_cog_orig);
  error_code |= p->deallocate_device(&gpu_buffers.d_ref_pos_cog);
  error_code |= p->deallocate_device(&gpu_buffers.d_com_cog_tbcount);
  error_code |= p->deallocate_device(&calc_fit_gradients_gpu_info.d_atom_grad);
  error_code |= p->deallocate_device(&calc_fit_gradients_gpu_info.d_sum_dxdq);
  error_code |= p->deallocate_device(&calc_fit_gradients_gpu_info.d_tbcount);
  error_code |= p->deallocate_device(&calc_fit_forces_gpu_info.d_atom_grad);
  error_code |= p->deallocate_device(&calc_fit_forces_gpu_info.d_sum_dxdq);
  error_code |= p->deallocate_device(&calc_fit_forces_gpu_info.d_tbcount);
  // error_code |= checkGPUError(cudaStreamDestroy(stream_ag_force));
  error_code |= p->deallocate_host(&h_sum_applied_colvar_force);
  error_code |= p->deallocate_host(&gpu_buffers.h_com);
  error_code |= p->deallocate_host(&gpu_buffers.h_cog);
  error_code |= p->deallocate_host(&gpu_buffers.h_cog_orig);
  if (debug_graphs.graph_calc_required_properties) {
    error_code |= checkGPUError(cudaGraphDestroy(
      debug_graphs.graph_calc_required_properties));
  }
  if (debug_graphs.graph_exec_calc_required_properties) {
    error_code |= checkGPUError(cudaGraphExecDestroy(
      debug_graphs.graph_exec_calc_required_properties));
  }
  debug_graphs.initialized = false;
  return error_code;
}

int cvm::atom_group::sync_to_gpu_buffers() {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvm::main()->proxy;
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_index, this->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_charge, this->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_mass, this->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_weight, this->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_pos, 3 * this->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_vel, 3 * this->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_total_force, 3 * this->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_grad, 3 * this->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_applied_force, 3 * this->num_atoms);
  error_code |= p->copy_HtoD(this->atoms_index.data(), this->gpu_buffers.d_atoms_index, this->num_atoms);
  error_code |= p->copy_HtoD(this->atoms_charge.data(), this->gpu_buffers.d_atoms_charge, this->num_atoms);
  error_code |= p->copy_HtoD(this->atoms_mass.data(), this->gpu_buffers.d_atoms_mass, this->num_atoms);
  error_code |= p->copy_HtoD(this->atoms_pos.data(), this->gpu_buffers.d_atoms_pos, 3 * this->num_atoms);
  error_code |= p->copy_HtoD(this->atoms_vel.data(), this->gpu_buffers.d_atoms_vel, 3 * this->num_atoms);
  error_code |= p->copy_HtoD(this->atoms_grad.data(), this->gpu_buffers.d_atoms_grad, 3 * this->num_atoms);
  error_code |= p->copy_HtoD(this->atoms_total_force.data(), this->gpu_buffers.d_atoms_total_force, 3 * this->num_atoms);
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
#endif
  return error_code;
}

int cvm::atom_group::clear_gpu_buffers() {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvm::main()->proxy;
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
  p->clear_device_array(gpu_buffers.d_atoms_index, num_atoms);
  p->clear_device_array(gpu_buffers.d_atoms_mass, num_atoms);
  p->clear_device_array(gpu_buffers.d_atoms_charge, num_atoms);
  p->clear_device_array(gpu_buffers.d_atoms_weight, num_atoms);
  p->clear_device_array(gpu_buffers.d_atoms_pos, 3 * num_atoms);
  p->clear_device_array(gpu_buffers.d_atoms_grad, 3 * num_atoms);
  p->clear_device_array(gpu_buffers.d_atoms_total_force, 3 * num_atoms);
  p->clear_device_array(gpu_buffers.d_atoms_vel, 3 * num_atoms);
  p->clear_device_array(gpu_buffers.d_atoms_applied_force, 3 * num_atoms);
#elif defined(COLVARS_SYCL)
// TODO: COLVARS_SYCL
#endif
  return error_code;
}

int cvm::atom_group::add_reset_atoms_data_nodes(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  if (!is_enabled(f_ag_scalable)) {
#define ADD_CLEAR_FIELD_NODE(fieldName, numAtoms) do {\
cudaGraphNode_t clear_ ## fieldName ## _node ;\
error_code |= colvars_gpu::add_clear_array_node( \
  gpu_buffers.d_ ## fieldName, 3 * numAtoms, \
  clear_ ## fieldName ## _node , graph, {});\
  nodes_map[COLVARS_STRINGIFY(clear_ ## fieldName)] = clear_ ## fieldName ## _node;\
} while (0);
    // ADD_CLEAR_FIELD_NODE(atoms_pos, num_atoms);
    // ADD_CLEAR_FIELD_NODE(atoms_vel, num_atoms);
    ADD_CLEAR_FIELD_NODE(atoms_grad, num_atoms);
    ADD_CLEAR_FIELD_NODE(atoms_total_force, num_atoms);
    ADD_CLEAR_FIELD_NODE(atoms_applied_force, num_atoms);
#undef ADD_CLEAR_FIELD_NODE
  }
  return error_code;
}

int cvm::atom_group::add_read_positions_nodes(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  if (b_dummy) return error_code;
  colvarproxy *p = cvm::main()->proxy;
  if (!is_enabled(f_ag_scalable)) {
    cudaGraphNode_t read_positions_node;
    std::vector<cudaGraphNode_t> dependencies;
    // We must wait for the clearing of device arrays
    // ADD_DEPENDENCY(clear_atoms_pos, dependencies, nodes_map);
    error_code |= colvars_gpu::atoms_pos_from_proxy(
      gpu_buffers.d_atoms_index, p->proxy_atoms_positions_gpu(),
      gpu_buffers.d_atoms_pos, num_atoms, p->get_atom_ids()->size(),
      read_positions_node,
      graph, dependencies);
    nodes_map["read_positions"] = read_positions_node;
  }
  if (fitting_group) {
    auto* current_fitting_group = fitting_group;
    // TODO: Do we really allow nesting fitting groups?
    int num_nested_fitting_group = 0;
    while (current_fitting_group) {
      cudaGraphNode_t read_fitting_group_positions_node;
      // TODO: Does the reading of the fitting group positions
      // depend on any other operations?
      error_code |= colvars_gpu::atoms_pos_from_proxy(
        current_fitting_group->gpu_buffers.d_atoms_index,
        p->proxy_atoms_positions_gpu(),
        current_fitting_group->gpu_buffers.d_atoms_pos,
        current_fitting_group->num_atoms,
        p->get_atom_ids()->size(),
        read_fitting_group_positions_node,
        graph, {});
      nodes_map["read_fitting_group_positions"] =
        read_fitting_group_positions_node;
      current_fitting_group = current_fitting_group->fitting_group;
      num_nested_fitting_group++;
    }
    if (num_nested_fitting_group > 1) {
      return cvm::error("Nesting fitting groups are not supported.\n");
    }
  }
  return error_code;
}

int cvm::atom_group::add_calc_required_properties_nodes(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
  const std::vector<cudaGraphNode_t>& extra_initial_dependencies) {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvm::main()->proxy;
  // COM and COG of the main group
  if (b_dummy) {
    com = dummy_atom_pos;
    cog = dummy_atom_pos;
    if (cvm::debug()) {
      cvm::log("Dummy atom center of mass = "+cvm::to_str(com)+"\n");
      cvm::log("Dummy atom center of cog = "+cvm::to_str(cog)+"\n");
    }
    // TODO: Do I need to copy dummy_atom_pos to GPU?
  } else if (is_enabled(f_ag_scalable)) {
    // com = (cvm::proxy)->get_atom_group_com(index);
    error_code |= cvm::error("BUG: GPU buffers are not implemented with scalable atom group.\n");
  } else {
    // Reset center-of-mass
    cudaGraphNode_t reset_com_node;
    error_code |= colvars_gpu::add_clear_array_node(
      gpu_buffers.d_com, 1, reset_com_node, graph, {});
    nodes_map["reset_com"] = reset_com_node;
    // Reset center-of-geometry
    cudaGraphNode_t reset_cog_node;
    error_code |= colvars_gpu::add_clear_array_node(
      gpu_buffers.d_cog, 1, reset_cog_node, graph, {});
    nodes_map["reset_cog"] = reset_cog_node;
    // Add kernel node
    cudaGraphNode_t calc_com_cog_node;
    std::vector<cudaGraphNode_t> dependencies = extra_initial_dependencies;
    ADD_DEPENDENCY_IF(read_positions, dependencies, nodes_map);
    // ADD_DEPENDENCY(reset_com_cog_tbcounter, dependencies, nodes_map);
    ADD_DEPENDENCY(reset_com, dependencies, nodes_map);
    ADD_DEPENDENCY(reset_cog, dependencies, nodes_map);
    error_code |= colvars_gpu::atoms_calc_cog_com(
      gpu_buffers.d_atoms_pos, gpu_buffers.d_atoms_mass, num_atoms,
      gpu_buffers.d_cog, gpu_buffers.d_com,
      gpu_buffers.h_cog, gpu_buffers.h_com, total_mass,
      gpu_buffers.d_com_cog_tbcount,
      calc_com_cog_node, graph, dependencies);
    nodes_map["calc_com_cog"] = calc_com_cog_node;
  }

  // Fitting group cog
  if (!is_enabled(f_ag_scalable)) {
    if (is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) {
      if (fitting_group) {
        if (fitting_group->b_dummy) {
          fitting_group->cog = fitting_group->dummy_atom_pos;
        } else {
          // Reset fitting group COG
          cudaGraphNode_t reset_fitting_group_cog_node;
          error_code |= colvars_gpu::add_clear_array_node(
            fitting_group->gpu_buffers.d_cog, 1,
            reset_fitting_group_cog_node,
            graph, {});
          nodes_map["reset_fitting_group_cog"] = reset_fitting_group_cog_node;
          // Add COG kernel node
          std::vector<cudaGraphNode_t> dependencies = extra_initial_dependencies;
          ADD_DEPENDENCY_IF(read_fitting_group_positions, dependencies, nodes_map);
          ADD_DEPENDENCY(reset_fitting_group_cog, dependencies, nodes_map);
          // ADD_DEPENDENCY(reset_fitting_group_cog_tbcounter, dependencies, nodes_map);
          cudaGraphNode_t calc_fitting_group_cog_node;
          error_code |= colvars_gpu::atoms_calc_cog(
            fitting_group->gpu_buffers.d_atoms_pos,
            fitting_group->num_atoms,
            fitting_group->gpu_buffers.d_cog,
            fitting_group->gpu_buffers.h_cog,
            fitting_group->gpu_buffers.d_com_cog_tbcount,
            calc_fitting_group_cog_node,
            graph,
            dependencies);
          nodes_map["calc_fitting_group_cog"] = calc_fitting_group_cog_node;
        }
      }

      // calc_apply_roto_translation()
      // store the laborarory-frame COGs for when they are needed later
      // Copy d_cog to d_cog_orig
      std::vector<cudaGraphNode_t> dependencies_cog_orig;
      ADD_DEPENDENCY(calc_com_cog, dependencies_cog_orig, nodes_map);
      cudaGraphNode_t save_cog_orig_node;
      error_code |= colvars_gpu::add_copy_node(
        gpu_buffers.d_cog, gpu_buffers.d_cog_orig, 1, cudaMemcpyDeviceToDevice,
        save_cog_orig_node,
        graph, dependencies_cog_orig);
      nodes_map["save_cog_orig"] = save_cog_orig_node;
      if (fitting_group) {
        std::vector<cudaGraphNode_t> dependencies_fitting_group_cog_orig;
        ADD_DEPENDENCY(calc_fitting_group_cog, dependencies_fitting_group_cog_orig, nodes_map);
        cudaGraphNode_t save_fitting_group_cog_orig_node;
        error_code |= colvars_gpu::add_copy_node(
          fitting_group->gpu_buffers.d_cog,
          fitting_group->gpu_buffers.d_cog_orig,
          1, cudaMemcpyDeviceToDevice,
          save_fitting_group_cog_orig_node,
          graph,
          dependencies_fitting_group_cog_orig);
        nodes_map["save_fitting_group_cog_orig"] = save_fitting_group_cog_orig_node;
      }

      // center on the origin first
      if (is_enabled(f_ag_center)) {
        const cvm::rvector* d_rpg_cog =
          fitting_group ? fitting_group->gpu_buffers.d_cog : this->gpu_buffers.d_cog;
        // Apply translation
        std::vector<cudaGraphNode_t> dependencies;
        ADD_DEPENDENCY(save_cog_orig, dependencies, nodes_map);
        cudaGraphNode_t move_to_origin_node;
        error_code |= colvars_gpu::apply_translation(
          gpu_buffers.d_atoms_pos, -1.0, d_rpg_cog,
          num_atoms, move_to_origin_node,
          graph, dependencies
        );
        nodes_map["move_to_origin"] = move_to_origin_node;
        if (fitting_group) {
          std::vector<cudaGraphNode_t> dependencies_fitting_group_translate;
          ADD_DEPENDENCY(save_fitting_group_cog_orig,
                          dependencies_fitting_group_translate,
                          nodes_map);
          cudaGraphNode_t move_fitting_to_origin_node;
          error_code |= colvars_gpu::apply_translation(
            fitting_group->gpu_buffers.d_atoms_pos, -1.0, d_rpg_cog,
            fitting_group->num_atoms,
            move_fitting_to_origin_node,
            graph,
            dependencies_fitting_group_translate);
          nodes_map["move_fitting_to_origin"] = move_fitting_to_origin_node;
        }
      }

      // Save the unrotated frame for fit gradients
      if (is_enabled(f_ag_fit_gradients) && !b_dummy) {
        std::vector<cudaGraphNode_t> dependencies;
        ADD_DEPENDENCY(move_to_origin, dependencies, nodes_map);
        cudaGraphNode_t copy_unrotated_positions_node;
        error_code |= colvars_gpu::add_copy_node(
          gpu_buffers.d_atoms_pos, gpu_buffers.d_atoms_pos_unrotated,
          3 * num_atoms, cudaMemcpyDeviceToDevice,
          copy_unrotated_positions_node,
          graph, dependencies);
        nodes_map["copy_unrotated_positions"] = copy_unrotated_positions_node;
      }

      // rotate the group (around the center of geometry if f_ag_center is
      // enabled, around the origin otherwise)
      if (is_enabled(f_ag_rotate)) {
        auto* group_for_fit = fitting_group ? fitting_group : this;
        error_code |= rot_gpu.add_optimal_rotation_nodes(
          group_for_fit->gpu_buffers.d_atoms_pos,
          gpu_buffers.d_ref_pos, group_for_fit->size(),
          num_ref_pos, graph, nodes_map);
        // Rotate atoms
        std::vector<cudaGraphNode_t> dependencies_rotate;
        ADD_DEPENDENCY(calc_optimal_rotation, dependencies_rotate, nodes_map);
        // If we have f_ag_fit_gradients, then we need to wait for the copying of unrotated positions
        ADD_DEPENDENCY_IF(copy_unrotated_positions, dependencies_rotate, nodes_map);
        // The atom group may be moved to origin (but not always)
        ADD_DEPENDENCY_IF(move_to_origin, dependencies_rotate, nodes_map);
        cudaGraphNode_t rotate_node;
        error_code |= colvars_gpu::rotate_with_quaternion(
          gpu_buffers.d_atoms_pos, rot_gpu.get_q(), num_atoms,
          rotate_node, graph,
          dependencies_rotate);
        nodes_map["rotate"] = rotate_node;
        if (fitting_group) {
          std::vector<cudaGraphNode_t> dependencies_rotate_fitting_group;
          ADD_DEPENDENCY(calc_optimal_rotation,
                         dependencies_rotate_fitting_group,
                         nodes_map);
          ADD_DEPENDENCY_IF(move_fitting_to_origin,
                            dependencies_rotate_fitting_group,
                            nodes_map);
          cudaGraphNode_t rotate_fitting_group_node;
          error_code |= colvars_gpu::rotate_with_quaternion(
            fitting_group->gpu_buffers.d_atoms_pos,
            rot_gpu.get_q(),
            fitting_group->num_atoms,
            rotate_fitting_group_node,
            graph,
            dependencies_rotate_fitting_group);
          nodes_map["rotate_fitting_group"] = rotate_fitting_group_node;
        }
      }

      // align with the center of geometry of ref_pos
      if (is_enabled(f_ag_center) && !is_enabled(f_ag_center_origin)) {
        std::vector<cudaGraphNode_t> dependencies_move_to_ref_cog;
        ADD_DEPENDENCY(calc_com_cog, dependencies_move_to_ref_cog, nodes_map);
        // TODO: It looks like that the moving to COG of ref_pos can depend on
        // any of the following operations, so I add them all. Is that right??
        ADD_DEPENDENCY_IF(rotate, dependencies_move_to_ref_cog, nodes_map);
        ADD_DEPENDENCY_IF(copy_unrotated_positions, dependencies_move_to_ref_cog, nodes_map);
        ADD_DEPENDENCY_IF(move_to_origin, dependencies_move_to_ref_cog, nodes_map);
        cudaGraphNode_t move_to_ref_cog_node;
        error_code |= colvars_gpu::apply_translation(
          gpu_buffers.d_atoms_pos, 1.0, gpu_buffers.d_ref_pos_cog, num_atoms,
          move_to_ref_cog_node, graph,
          dependencies_move_to_ref_cog);
        nodes_map["move_to_ref_cog"] = move_to_ref_cog_node;
        if (fitting_group) {
          std::vector<cudaGraphNode_t> dependencies_move_fitting_group_to_ref_cog;
          ADD_DEPENDENCY_IF(
            calc_fitting_group_cog,
            dependencies_move_fitting_group_to_ref_cog,
            nodes_map);
          ADD_DEPENDENCY_IF(
            move_fitting_to_origin,
            dependencies_move_fitting_group_to_ref_cog,
            nodes_map);
          ADD_DEPENDENCY_IF(
            rotate_fitting_group,
            dependencies_move_fitting_group_to_ref_cog,
            nodes_map);
          cudaGraphNode_t move_fitting_group_to_ref_cog_node;
          error_code |= colvars_gpu::apply_translation(
            fitting_group->gpu_buffers.d_atoms_pos, 1.0,
            gpu_buffers.d_ref_pos_cog, num_atoms,
            move_fitting_group_to_ref_cog_node,
            graph,
            dependencies_move_fitting_group_to_ref_cog);
          nodes_map["move_fitting_group_to_ref_cog"] =
            move_fitting_group_to_ref_cog_node;
        }
      }

      // update COM and COG after fitting
      std::vector<cudaGraphNode_t> dependencies;
      // We reuse the atomic thread block counter so we need to wait for calc_com_cog
      ADD_DEPENDENCY(calc_com_cog, dependencies, nodes_map);
      // Reset center-of-mass
      // COM or COG may be used in any of the following operations:
      //   move_to_origin, move_fitting_to_origin, save_cog_orig
      // TODO: Check if I still miss something?
      ADD_DEPENDENCY_IF(move_to_origin, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(move_fitting_to_origin, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(save_cog_orig, dependencies, nodes_map);
      cudaGraphNode_t reset_com2_node;
      error_code |= colvars_gpu::add_clear_array_node(
        gpu_buffers.d_com, 1, reset_com2_node, graph, dependencies);
      nodes_map["reset_com2"] = reset_com2_node;
      // Reset center-of-geometry
      cudaGraphNode_t reset_cog2_node;
      error_code |= colvars_gpu::add_clear_array_node(
        gpu_buffers.d_cog, 1, reset_cog2_node, graph, dependencies);
      nodes_map["reset_cog2"] = reset_cog2_node;
      // Re-calculate COM and COG
      // ADD_DEPENDENCY(reset_com_cog_tbcounter2, dependencies, nodes_map);
      ADD_DEPENDENCY(reset_com2, dependencies, nodes_map);
      ADD_DEPENDENCY(reset_cog2, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(rotate, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(move_to_ref_cog, dependencies, nodes_map);
      cudaGraphNode_t calc_com_cog2_node;
      error_code |= colvars_gpu::atoms_calc_cog_com(
        gpu_buffers.d_atoms_pos, gpu_buffers.d_atoms_mass, num_atoms,
        gpu_buffers.d_cog, gpu_buffers.d_com,
        gpu_buffers.h_cog, gpu_buffers.h_com,
        total_mass, gpu_buffers.d_com_cog_tbcount,
        calc_com_cog2_node, graph, dependencies);
      nodes_map["calc_com_cog2"] = calc_com_cog2_node;
      if (fitting_group) {
        dependencies.clear();
        // fitting_group->d_com_cog_tbcount is used in calc_fitting_group_cog
        // so we have to wait for the previous operation.
        ADD_DEPENDENCY(calc_fitting_group_cog, dependencies, nodes_map);
        // The following operations may or may not exist,
        // but if any of them exist, we need to wait for them before
        // resetting the fitting group COG.
        ADD_DEPENDENCY_IF(save_fitting_group_cog_orig, dependencies, nodes_map);
        ADD_DEPENDENCY_IF(move_fitting_to_origin, dependencies, nodes_map);
        cudaGraphNode_t reset_fitting_group_cog2_node;
        error_code |= colvars_gpu::add_clear_array_node(
          fitting_group->gpu_buffers.d_cog, 1,
          reset_fitting_group_cog2_node,
          graph, dependencies);
        nodes_map["reset_fitting_group_cog2"] = reset_fitting_group_cog2_node;
        // Re-calculate fitting group COG
        // ADD_DEPENDENCY(
        //   reset_fitting_group_cog_tbcounter2, dependencies, nodes_map);
        ADD_DEPENDENCY(reset_fitting_group_cog2, dependencies, nodes_map);
        ADD_DEPENDENCY_IF(rotate_fitting_group, dependencies, nodes_map);
        ADD_DEPENDENCY_IF(move_fitting_group_to_ref_cog, dependencies, nodes_map);
        cudaGraphNode_t calc_fitting_group_cog2_node;
        error_code |= colvars_gpu::atoms_calc_cog(
          fitting_group->gpu_buffers.d_atoms_pos,
          fitting_group->num_atoms,
          fitting_group->gpu_buffers.d_cog,
          fitting_group->gpu_buffers.h_cog,
          fitting_group->gpu_buffers.d_com_cog_tbcount,
          calc_fitting_group_cog2_node,
          graph, dependencies);
        nodes_map["calc_fitting_group_cog2"] = calc_fitting_group_cog2_node;
      }
    }
  }

  return error_code;
}

int cvm::atom_group::add_update_cpu_buffers_nodes(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  if (b_dummy) return error_code;
  colvarproxy *p = cvm::main()->proxy;
  if (!is_enabled(f_ag_scalable)) {
    std::vector<cudaGraphNode_t> dependencies;
    ADD_DEPENDENCY_IF(read_positions, dependencies, nodes_map);
    ADD_DEPENDENCY_IF(move_to_origin, dependencies, nodes_map);
    ADD_DEPENDENCY_IF(rotate, dependencies, nodes_map);
    ADD_DEPENDENCY_IF(move_to_ref_cog, dependencies, nodes_map);
    cudaGraphNode_t copy_atoms_to_host_node;
    error_code |= colvars_gpu::add_copy_node(
      gpu_buffers.d_atoms_pos, atoms_pos.data(), 3 * num_atoms,
      cudaMemcpyDeviceToHost, copy_atoms_to_host_node,
      graph, dependencies);
    nodes_map["copy_atoms_to_host"] = copy_atoms_to_host_node;
    if (fitting_group) {
      dependencies.clear();
      ADD_DEPENDENCY_IF(read_fitting_group_positions, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(move_fitting_to_origin, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(rotate_fitting_group, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(move_fitting_group_to_ref_cog, dependencies, nodes_map);
      cudaGraphNode_t copy_fitting_group_atoms_to_host_node;
      error_code |= colvars_gpu::add_copy_node(
        fitting_group->gpu_buffers.d_atoms_pos,
        fitting_group->atoms_pos.data(), 3 * num_atoms,
        cudaMemcpyDeviceToHost,
        copy_fitting_group_atoms_to_host_node,
        graph, dependencies);
      nodes_map["copy_fitting_group_atoms_to_host"] = copy_fitting_group_atoms_to_host_node;
    }
    if (is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) {
      if (is_enabled(f_ag_fit_gradients) && !b_dummy) {
        atoms_pos_unrotated.resize(3 * num_atoms);
        dependencies.clear();
        ADD_DEPENDENCY_IF(copy_unrotated_positions, dependencies, nodes_map);
        cudaGraphNode_t copy_unrotated_positions_to_host_node;
        error_code |= colvars_gpu::add_copy_node(
          gpu_buffers.d_atoms_pos_unrotated,
          atoms_pos_unrotated.data(), 3 * num_atoms,
          cudaMemcpyDeviceToHost,
          copy_unrotated_positions_to_host_node,
          graph, dependencies);
        nodes_map["copy_unrotated_positions_to_host"] = copy_unrotated_positions_to_host_node;
      }
    }
  }
  return error_code;
}

int cvm::atom_group::after_read_data_sync(
  bool copy_to_cpu, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  // Update the COM
  if (b_dummy) {
    com = dummy_atom_pos;
    if (cvm::debug()) {
      cvm::log("Dummy atom center of mass = "+cvm::to_str(com)+"\n");
    }
  } else if (is_enabled(f_ag_scalable)) {
    com = (cvm::proxy)->get_atom_group_com(index);
  } else {
    com.reset();
    if (copy_to_cpu) {
      com = *(gpu_buffers.h_com);
    }
  }
  // Update the COG
  if (b_dummy) {
    cog = dummy_atom_pos;
  } else {
    cog.reset();
    if (copy_to_cpu) {
      cog = *(gpu_buffers.h_cog);
    }
  }

  if (!is_enabled(f_ag_scalable)) {
    if (is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) {
      if (fitting_group) {
        // Update the fitting group COG
        if (fitting_group->b_dummy) {
          fitting_group->cog = fitting_group->dummy_atom_pos;
        } else {
          fitting_group->cog.reset();
          if (copy_to_cpu) {
            fitting_group->cog = *(fitting_group->gpu_buffers.h_cog);
          }
        }
      }

      if (copy_to_cpu) {
        cog_orig = *(gpu_buffers.h_cog_orig);
      }
      if (fitting_group) {
        if (copy_to_cpu) {
          fitting_group->cog_orig = *(fitting_group->gpu_buffers.h_cog_orig);
        }
      }

      if (is_enabled(f_ag_rotate)) {
        if (copy_to_cpu) rot_gpu.to_cpu(rot);
        rot_gpu.after_sync_check();
      }
    }
  }
  return error_code;
}

int cvm::atom_group::add_calc_fit_gradients_nodes(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
  bool use_cpu_buffers) {
  int error_code = COLVARS_OK;
  if (b_dummy || ! is_enabled(f_ag_fit_gradients)) return error_code;
  cvm::atom_group *group_for_fit = fitting_group ? fitting_group : this;
  colvarproxy* p = cvm::main()->proxy;
  // First, clear the temporary variables
  cudaGraphNode_t clear_atoms_grad_node;
  cudaGraphNode_t clear_sum_dxdq_node;
  // cudaGraphNode_t clear_tbcount_node;
  error_code |= colvars_gpu::add_clear_array_node(
    calc_fit_gradients_gpu_info.d_atom_grad,
    1, clear_atoms_grad_node, graph, {});
  error_code |= colvars_gpu::add_clear_array_node(
    calc_fit_gradients_gpu_info.d_sum_dxdq,
    1, clear_sum_dxdq_node, graph, {});
  // error_code |= colvars_gpu::add_clear_array_node(
  //   calc_fit_gradients_gpu_info.d_tbcount,
  //   1, clear_tbcount_node, graph, {});
  nodes_map["clear_atoms_grad"] = clear_atoms_grad_node;
  nodes_map["clear_sum_dxdq"] = clear_sum_dxdq_node;
  // nodes_map["clear_tbcount"] = clear_tbcount_node;
  // If the CVC updates the gradients on CPU, then we need to copy them to GPU
  if (use_cpu_buffers) {
    cudaGraphNode_t copy_grad_HtoD_node;
    error_code |= colvars_gpu::add_copy_node(
      atoms_grad.data(), gpu_buffers.d_atoms_grad, 3 * num_atoms,
      cudaMemcpyHostToDevice, copy_grad_HtoD_node, graph,
      {clear_atoms_grad_node});
    nodes_map["copy_grad_HtoD"] = copy_grad_HtoD_node;
  }
  // Prepare the rotation derivative at the same time
  if (is_enabled(f_ag_rotate) && rot_deriv_gpu) {
    error_code |= rot_deriv_gpu->add_prepare_derivative_nodes(
      rotation_derivative_dldq::use_dq, graph, nodes_map);
  }
  // Loop over the main group atoms (does not require rot_deriv_gpu)
  if (is_enabled(f_ag_rotate) || is_enabled(f_ag_center)) {
    std::vector<cudaGraphNode_t> dependencies_main;
    ADD_DEPENDENCY(clear_atoms_grad, dependencies_main, nodes_map);
    ADD_DEPENDENCY(clear_sum_dxdq, dependencies_main, nodes_map);
    // ADD_DEPENDENCY(clear_tbcount, dependencies_main, nodes_map);
    ADD_DEPENDENCY_IF(copy_grad_HtoD, dependencies_main, nodes_map);
    cudaGraphNode_t calc_fit_forces_loop1_node;
    error_code |= colvars_gpu::calc_fit_gradients_impl_loop1(
      gpu_buffers.d_atoms_pos_unrotated, gpu_buffers.d_atoms_grad, rot_gpu.get_q(),
      num_atoms, group_for_fit->size(),
      calc_fit_gradients_gpu_info.d_atom_grad,
      calc_fit_gradients_gpu_info.d_sum_dxdq,
      calc_fit_gradients_gpu_info.d_tbcount,
      is_enabled(f_ag_center),
      is_enabled(f_ag_rotate),
      calc_fit_forces_loop1_node,
      graph, dependencies_main);
    nodes_map["calc_fit_gradients_loop1"] = calc_fit_forces_loop1_node;
    // Loop over the fitting group
    // if (gpu_buffers.d_fit_gradients == nullptr) {
    //   // TODO: Avoid allocation here
    //   error_code |= p->allocate_device(&gpu_buffers.d_fit_gradients, 3 * group_for_fit->size());
    //   gpu_buffers.d_fit_gradients_size = group_for_fit->size();
    // }
    cudaGraphNode_t calc_fit_forces_loop2_node;
    std::vector<cudaGraphNode_t> dependencies_fit_gradients;
    ADD_DEPENDENCY(calc_fit_gradients_loop1, dependencies_fit_gradients, nodes_map);
    ADD_DEPENDENCY_IF(prepare_rotation_derivative, dependencies_fit_gradients, nodes_map);
    error_code |= colvars_gpu::calc_fit_gradients_impl_loop2(
      group_for_fit->gpu_buffers.d_fit_gradients, rot_deriv_gpu,
      calc_fit_gradients_gpu_info.d_atom_grad,
      calc_fit_gradients_gpu_info.d_sum_dxdq,
      group_for_fit->size(),
      is_enabled(f_ag_center),
      is_enabled(f_ag_rotate),
      calc_fit_forces_loop2_node,
      graph, dependencies_fit_gradients);
    nodes_map["calc_fit_gradients_loop2"] = calc_fit_forces_loop2_node;
    if (use_cpu_buffers) {
      if (group_for_fit->fit_gradients.empty()) {
        group_for_fit->fit_gradients.resize(3 * group_for_fit->size());
      }
      cudaGraphNode_t copy_fit_gradients_DtoH_node;
      std::vector<cudaGraphNode_t> dependencies_copy_fit_gradients;
      ADD_DEPENDENCY(calc_fit_gradients_loop2, dependencies_copy_fit_gradients, nodes_map);
      error_code |= colvars_gpu::add_copy_node(
        group_for_fit->gpu_buffers.d_fit_gradients, group_for_fit->fit_gradients.data(), 3 * group_for_fit->size(),
        cudaMemcpyDeviceToHost, copy_fit_gradients_DtoH_node, graph,
        dependencies_copy_fit_gradients);
      nodes_map["copy_fit_gradients_DtoH"] = copy_fit_gradients_DtoH_node;
    }
  }
  return error_code;
}

int cvm::atom_group::begin_apply_force_gpu() {
  int error_code = COLVARS_OK;
  h_sum_applied_colvar_force[0] = 0;
  use_apply_colvar_force = false;
  use_group_force = false;
  return error_code;
}

int cvm::atom_group::add_apply_force_nodes(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
  const std::vector<cudaGraphNode_t>& extra_initial_dependencies) {
  int error_code = COLVARS_OK;
  if (b_dummy) return error_code;
  colvarproxy* p = cvm::main()->proxy;
  // Check if any of the parent CVCs require CPU buffers
  // Get all parents
  const std::vector<colvardeps*> parents = get_parents();
  std::vector<colvar::cvc*> parent_components;
  // Cast the parents to colvarcomp objects if possible
  for (size_t i = 0; i < parents.size(); ++i) {
    try {
      parent_components.push_back(dynamic_cast<colvar::cvc*>(parents[i]));
    } catch (const std::bad_cast& exception) {
      // Ignore the bad cast
    }
  }
  auto check_cvc_cpu_buffers = [](const colvar::cvc* obj){return obj->is_enabled(f_cvc_require_cpu_buffers);};
  const bool any_require_cpu_buffers = std::any_of(parent_components.begin(), parent_components.end(), check_cvc_cpu_buffers);
  const bool all_require_cpu_buffers = std::all_of(parent_components.begin(), parent_components.end(), check_cvc_cpu_buffers);
  if (use_apply_colvar_force) {
    cvm::quaternion* q = nullptr;
    if (rot_gpu.initialized()) {
      q = rot_gpu.get_q();
    }
    if (any_require_cpu_buffers) {
      cudaGraphNode_t copy_grad_HtoD_node;
      error_code |= colvars_gpu::add_copy_node(
        atoms_grad.data(), gpu_buffers.d_atoms_grad, 3 * num_atoms,
        cudaMemcpyHostToDevice, copy_grad_HtoD_node, graph, {});
      nodes_map["copy_grad_HtoD"] = copy_grad_HtoD_node;
      if (!all_require_cpu_buffers) {
        const std::string error = "BUG: either none of the CVCs or "
          " all of the CVCs require CPU buffers!\n";
        return cvm::error(error);
      }
    }
    cudaGraphNode_t apply_colvar_force_to_proxy_node;
    std::vector<cudaGraphNode_t> dependencies = extra_initial_dependencies;
    ADD_DEPENDENCY_IF(copy_grad_HtoD, dependencies, nodes_map);
    error_code |= colvars_gpu::apply_colvar_force_to_proxy(
      gpu_buffers.d_atoms_index,
      p->proxy_atoms_new_colvar_forces_gpu(),
      gpu_buffers.d_atoms_grad,
      h_sum_applied_colvar_force,
      is_enabled(f_ag_rotate),
      q, num_atoms, p->get_atom_ids()->size(),
      apply_colvar_force_to_proxy_node, graph, dependencies);
    nodes_map["apply_colvar_force_to_proxy"] = apply_colvar_force_to_proxy_node;
  }
  if (use_group_force) {
    std::vector<cudaGraphNode_t> dependencies = extra_initial_dependencies;
    if (all_require_cpu_buffers) {
      // All CVCs using this group require CPU buffers, so
      // we can directly copy the data from group_forces to
      // d_atoms_applied_force.
      cudaGraphNode_t copy_forces_HtoD_node;
      error_code |= colvars_gpu::add_copy_node(
        group_forces.data(), gpu_buffers.d_atoms_applied_force, 3 * num_atoms,
        cudaMemcpyHostToDevice, copy_forces_HtoD_node, graph, dependencies);
      nodes_map["copy_forces_HtoD"] = copy_forces_HtoD_node;
    } else {
      if (any_require_cpu_buffers) {
        // Only some of the CVCs require the CPU buffers, so we
        // need to add the CPU buffers to d_atoms_applied_force
        cudaGraphNode_t accumulate_cpu_force_node;
        error_code |= colvars_gpu::accumulate_cpu_force(
          group_forces.data(), gpu_buffers.d_atoms_applied_force, num_atoms,
          accumulate_cpu_force_node, graph, dependencies);
        nodes_map["accumulate_cpu_force"] = accumulate_cpu_force_node;
      }
    }
    ADD_DEPENDENCY_IF(copy_forces_HtoD, dependencies, nodes_map);
    ADD_DEPENDENCY_IF(accumulate_cpu_force, dependencies, nodes_map);
    if (is_enabled(f_ag_rotate)) {
      // Rotate the forces back and add them to proxy
      cudaGraphNode_t apply_force_with_inverse_rotation_node;
      error_code |= colvars_gpu::apply_force_with_inverse_rotation(
        group_forces.data(), rot_gpu.get_q(), gpu_buffers.d_atoms_index,
        p->proxy_atoms_new_colvar_forces_gpu(), num_atoms,
        p->get_atom_ids()->size(),
        apply_force_with_inverse_rotation_node,
        graph, dependencies);
      nodes_map["apply_force_with_inverse_rotation"] =
        apply_force_with_inverse_rotation_node;
    } else {
      // Just add the forces to proxy
      cudaGraphNode_t apply_force_node;
      error_code |= colvars_gpu::apply_force(
        group_forces.data(), gpu_buffers.d_atoms_index,
        p->proxy_atoms_new_colvar_forces_gpu(), num_atoms,
        p->get_atom_ids()->size(), apply_force_node,
        graph, dependencies);
      nodes_map["apply_force"] = apply_force_node;
    }
    // ADD_DEPENDENCY_IF(apply_force_with_inverse_rotation, dependencies, nodes_map);
    // ADD_DEPENDENCY_IF(apply_force, dependencies, nodes_map);
    // dependencies = extra_initial_dependencies;
    if (is_enabled(f_ag_fit_gradients)) {
      auto* group_for_fit = this->fitting_group ? this->fitting_group : this;
      // Compute the forces on the fitting group and add them to proxy
      // Clear the temporary variables
      cudaGraphNode_t clear_atoms_grad_node;
      cudaGraphNode_t clear_sum_dxdq_node;
      cudaGraphNode_t clear_tbcount_node;
      error_code |= colvars_gpu::add_clear_array_node(
        calc_fit_forces_gpu_info.d_atom_grad,
        1, clear_atoms_grad_node, graph, {});
      error_code |= colvars_gpu::add_clear_array_node(
        calc_fit_forces_gpu_info.d_sum_dxdq,
        1, clear_sum_dxdq_node, graph, {});
      error_code |= colvars_gpu::add_clear_array_node(
        calc_fit_forces_gpu_info.d_tbcount,
        1, clear_tbcount_node, graph, {});
      nodes_map["clear_atoms_grad"] = clear_atoms_grad_node;
      nodes_map["clear_sum_dxdq"] = clear_sum_dxdq_node;
      nodes_map["clear_tbcount"] = clear_tbcount_node;
      ADD_DEPENDENCY_IF(clear_atoms_grad, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(clear_sum_dxdq, dependencies, nodes_map);
      ADD_DEPENDENCY_IF(clear_tbcount, dependencies, nodes_map);
      cudaGraphNode_t calc_fit_forces_loop1_node;
      error_code |= colvars_gpu::calc_fit_forces_impl_loop1(
        gpu_buffers.d_atoms_pos_unrotated,
        gpu_buffers.d_atoms_applied_force,
        rot_gpu.get_q(),
        num_atoms, group_for_fit->size(),
        calc_fit_forces_gpu_info.d_atom_grad,
        calc_fit_forces_gpu_info.d_sum_dxdq,
        calc_fit_forces_gpu_info.d_tbcount,
        is_enabled(f_ag_center),
        is_enabled(f_ag_rotate),
        calc_fit_forces_loop1_node,
        graph, dependencies);
      nodes_map["calc_fit_forces_loop1"] = calc_fit_forces_loop1_node;
      // Compute the forces on the fitting group and add them back to proxy
      cudaGraphNode_t calc_fit_forces_loop2_node;
      error_code |= colvars_gpu::calc_fit_forces_impl_loop2(
        rot_deriv_gpu,
        calc_fit_forces_gpu_info.d_atom_grad,
        calc_fit_forces_gpu_info.d_sum_dxdq,
        gpu_buffers.d_atoms_index, p->proxy_atoms_new_colvar_forces_gpu(),
        group_for_fit->size(), p->get_atom_ids()->size(),
        is_enabled(f_ag_center),
        is_enabled(f_ag_rotate),
        calc_fit_forces_loop2_node, graph,
        {calc_fit_forces_loop1_node});
      nodes_map["calc_fit_forces_loop2"] = calc_fit_forces_loop2_node;
    }
  }
  return error_code;
}

int cvm::atom_group::read_positions_gpu_debug(
  size_t change_atom_i, int xyz, bool to_cpu, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  colvarproxy *p = cvm::main()->proxy;
  error_code |= colvars_gpu::atoms_pos_from_proxy(
    gpu_buffers.d_atoms_index, p->proxy_atoms_positions_gpu(),
    gpu_buffers.d_atoms_pos, num_atoms, p->get_atom_ids()->size(),
    stream);
  error_code |= colvars_gpu::change_one_coordinate(
    gpu_buffers.d_atoms_pos, change_atom_i, xyz,
    cvm::debug_gradients_step_size, num_atoms, stream);
  if (to_cpu) {
    error_code |= p->copy_DtoH(
      gpu_buffers.d_atoms_pos, atoms_pos.data(), 3 * num_atoms);
    error_code |= checkGPUError(cudaStreamSynchronize(stream));
  }
  // error_code |= checkGPUError(cudaStreamSynchronize(stream));
  return error_code;
}

int cvm::atom_group::calc_required_properties_gpu_debug(bool to_cpu, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  if (!debug_graphs.initialized) {
    // Create the debug graph
    error_code |= checkGPUError(cudaGraphCreate(&debug_graphs.graph_calc_required_properties, 0));
    std::unordered_map<std::string, cudaGraphNode_t> nodes_map;
    error_code |= add_calc_required_properties_nodes(
      debug_graphs.graph_calc_required_properties, nodes_map);
    if (to_cpu) {
      error_code |= add_update_cpu_buffers_nodes(
        debug_graphs.graph_calc_required_properties, nodes_map);
    }
    error_code |= checkGPUError(cudaGraphInstantiate(
      &debug_graphs.graph_exec_calc_required_properties, debug_graphs.graph_calc_required_properties));
    debug_graphs.initialized = true;
  }
  error_code |= checkGPUError(cudaGraphLaunch(
    debug_graphs.graph_exec_calc_required_properties, stream));
  return error_code;
}

void cvm::atom_group::do_feature_side_effects_gpu(int id) {
  if (cvm::debug()) {
    cvm::log("cvm::atom_group::do_feature_side_effects_gpu.\n");
  }
  switch (id) {
    case f_ag_fit_gradients: {
      colvarproxy* p = cvm::main()->proxy;
      if (gpu_buffers.d_atoms_pos_unrotated == nullptr) {
        p->allocate_device(&gpu_buffers.d_atoms_pos_unrotated, 3 * num_atoms);
        gpu_buffers.d_atoms_pos_unrotated_size = num_atoms;
      }
      if (is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) {
        atom_group *group_for_fit = fitting_group ? fitting_group : this;
        if (group_for_fit->gpu_buffers.d_fit_gradients == nullptr) {
          p->allocate_device(&group_for_fit->gpu_buffers.d_fit_gradients,
                              3 * group_for_fit->size());
          p->clear_device_array(group_for_fit->gpu_buffers.d_fit_gradients, 3 * group_for_fit->size());
          if (cvm::debug()) {
            cvm::log("allocate d_fit_gradients at " + cvm::to_str((void*)group_for_fit->gpu_buffers.d_fit_gradients) + " size " + cvm::to_str(3 * group_for_fit->size()) + "\n");
          }
          group_for_fit->gpu_buffers.d_fit_gradients_size = group_for_fit->size();
        }
      }
      break;
    }
    case f_ag_rotate: {
      rot_gpu.init();
      break;
    }
  }
}
#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
