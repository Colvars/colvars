#ifndef COLVAR_GPU_CALC_H
#define COLVAR_GPU_CALC_H

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)

#include <vector>
#include "colvar_gpu_support.h"

class colvardeps;
class colvar;
class colvarmodule;

namespace colvars_gpu {
class colvarmodule_gpu_calc {
public:
  struct compute_node_t {
    // XXX TODO: There should be an independent data structure for AST
    // For the time being I have to hack the colvardeps which serves
    // partially as an AST (without any type information). There are two
    // kinds of dependencies, namely (i) the dependencies between computational
    // operations (cv-cvc, cvc-cvc and cvc-atom groups), and (ii) the
    // dependencies of different features. (i) should be similar to an AST
    // in a compiler, while (ii) should be a checker walking over the AST.
    colvardeps* colvar_node;
    cudaGraphNode_t child_graph_node;
    bool require_cpu_buffers;
  };
  class compute_gpu_graph_t {
  public:
    compute_gpu_graph_t();
    int init();
    ~compute_gpu_graph_t();
    void dump_graph(const char* filename);
    bool graph_exec_initialized;
    cudaGraph_t graph;
    cudaGraphExec_t graph_exec;
    std::vector<compute_node_t> nodes;
  };
  colvarmodule_gpu_calc();
  ~colvarmodule_gpu_calc() {}
  int init();
  int calc_cvs(const std::vector<colvar*>& colvars, colvarmodule* colvar_module);
  int apply_forces(const std::vector<colvar*>& colvars, colvarmodule* colvar_module);
private:
  compute_gpu_graph_t read_data_compute;
  compute_gpu_graph_t calc_fit_gradients_compute;
  compute_gpu_graph_t apply_forces_compute;
  /**
   * @note I have no way to forward declare cvm::atom_group here so I
   * have to use colvardeps as a workaround.
   */
  std::vector<colvardeps*> forced_atom_groups;
  // Timers for profiling
  #if defined (COLVARS_NVTX_PROFILING)
  colvars_gpu::colvar_nvtx_prof ag_read_data_prof;
  colvars_gpu::colvar_nvtx_prof cvc_calc_value_prof;
  colvars_gpu::colvar_nvtx_prof cvc_calc_gradients_prof;
  colvars_gpu::colvar_nvtx_prof ag_calc_fit_gradients_prof;
  colvars_gpu::colvar_nvtx_prof cvc_calc_Jacobian_derivative_prof;
  colvars_gpu::colvar_nvtx_prof cvc_calc_total_force_prof;
  colvars_gpu::colvar_nvtx_prof cv_collect_cvc_data_prof;
  colvars_gpu::colvar_nvtx_prof apply_forces_prof;
  #endif

  int cv_update_flags(const std::vector<colvar*>& colvars);
  int cvc_calc_total_force(
    const std::vector<colvar*>& colvars,
    colvarmodule* colvar_module,
    bool use_current_step = false);
  int atom_group_read_data_gpu(
    const std::vector<colvar*>& colvars,
    colvarmodule_gpu_calc::compute_gpu_graph_t& g,
    colvarmodule* colvar_module);
  int cvc_calc_value(
    const std::vector<colvar*>& colvars,
    colvarmodule* colvar_module);
  int cvc_calc_gradients(
    const std::vector<colvar*>& colvars,
    colvarmodule* colvar_module);
  int atom_group_calc_fit_gradients(
    const std::vector<colvar*>& colvars,
    colvarmodule_gpu_calc::compute_gpu_graph_t& g,
    colvarmodule* colvar_module);
  int cvc_debug_gradients(
    const std::vector<colvar*>& colvars,
    colvarmodule* colvar_module);
  int cvc_calc_Jacobian_derivative(
    const std::vector<colvar*>& colvars,
    colvarmodule* colvar_module);
  int cv_collect_cvc_data(
    const std::vector<colvar*>& colvars,
    colvarmodule* colvar_module);
};
}

#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
#endif // COLVAR_GPU_CALC_H
