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
  colvarmodule_gpu_calc() {}
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
};
}

#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
#endif // COLVAR_GPU_CALC_H
