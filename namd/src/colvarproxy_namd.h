// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_NAMD_H
#define COLVARPROXY_NAMD_H

#ifndef NAMD_VERSION_NUMBER
// Assume 2.14b1 for now until the NAMD macro is merged
#define NAMD_VERSION_NUMBER 34471681
#endif

#include <memory>

#include "colvarproxy_namd_version.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"


class Controller;
class GlobalMasterColvars;
class GridforceFullMainGrid;
class Random;
class SimParameters;


/// Communication between colvars and NAMD (implementation of \link colvarproxy \endlink)
class colvarproxy_namd : public colvarproxy {

protected:

  /// Pointer to the parent GlobalMaster object
  GlobalMasterColvars *globalmaster = nullptr;

  /// \brief Array of atom indices (relative to the colvarproxy arrays),
  /// usedfor faster copy of atomic data
  std::vector<int> atoms_map;

  /// Pointer to the NAMD simulation input object
  SimParameters *simparams;

  /// Pointer to Controller object
  Controller const *controller;

  /// NAMD-style PRNG object
  std::unique_ptr<Random> random;

  bool first_timestep;
  cvm::step_number previous_NAMD_step;

  /// Used to submit restraint energy as MISC
#if !defined (NAMD_UNIFIED_REDUCTION)
  SubmitReduction *reduction;
#endif
#if defined(NODEGROUP_FORCE_REGISTER) && !defined(NAMD_UNIFIED_REDUCTION)
  NodeReduction *nodeReduction;
#endif

  /// Accelerated MD reweighting factor
  bool accelMDOn;
  cvm::real amd_weight_factor;
  void update_accelMD_info();

public:

  void init_tcl_pointers() override;

  colvarproxy_namd(GlobalMasterColvars *gm);
  ~colvarproxy_namd();

  int setup() override;
  int reset() override;

  /// Get the target temperature from the NAMD thermostats supported so far
  int update_target_temperature();

  /// Allocate an atoms map with the same size as the NAMD topology
  void init_atoms_map();

  // synchronize the local arrays with requested or forced atoms
  int update_atoms_map(AtomIDList::const_iterator begin,
                       AtomIDList::const_iterator end);

  void calculate();

  void log(std::string const &message) override;
  void error(std::string const &message) override;
  int set_unit_system(std::string const &units_in, bool check_only) override;
  void add_energy(cvm::real energy) override;
  void request_total_force(bool yesno) override;

  bool total_forces_enabled() const override
  {
    return total_force_requested;
  }

  int run_force_callback() override;
  int run_colvar_callback(std::string const &name,
                          std::vector<const colvarvalue *> const &cvcs,
                          colvarvalue &value) override;
  int run_colvar_gradient_callback(std::string const &name,
                                   std::vector<const colvarvalue *> const &cvcs,
                                   std::vector<cvm::matrix2d<cvm::real> > &gradient) override;

  cvm::real rand_gaussian() override;

  cvm::real get_accelMD_factor() const override;

  bool accelMD_enabled() const override;

#if CMK_SMP && USE_CKLOOP
  colvarproxy::smp_mode_t get_smp_mode() const override;

  int set_smp_mode(smp_mode_t mode) override;

  smp_mode_t get_preferred_smp_mode() const override {
    return smp_mode_t::cvcs;
  }
  std::vector<smp_mode_t> get_available_smp_modes() const override {
    std::vector<colvarproxy_smp::smp_mode_t> available_modes{
      smp_mode_t::cvcs,
      smp_mode_t::inner_loop,
      smp_mode_t::none
    };
    return available_modes;
  }

  int smp_loop(int n_items, std::function<int (int)> const &worker) override;

  int smp_biases_loop() override;

  int smp_biases_script_loop() override;

  friend void calc_colvars_items_smp(int first, int last, void *result, int paramNum, void *param);
  friend void calc_cv_biases_smp(int first, int last, void *result, int paramNum, void *param);
  friend void calc_cv_scripted_forces(int paramNum, void *param);

  int smp_thread_id()
  {
    return CkMyRank();
  }

  int smp_num_threads()
  {
    return CkMyNodeSize();
  }

protected:

  CmiNodeLock charm_lock_state;

public:

  int smp_lock()
  {
    CmiLock(charm_lock_state);
    return COLVARS_OK;
  }

  int smp_trylock()
  {
    const int ret = CmiTryLock(charm_lock_state);
    if (ret == 0) return COLVARS_OK;
    else return COLVARS_ERROR;
  }

  int smp_unlock()
  {
    CmiUnlock(charm_lock_state);
    return COLVARS_OK;
  }

#endif // #if CMK_SMP && USE_CKLOOP

  int check_replicas_enabled() override;
  int replica_index() override;
  int num_replicas() override;
  void replica_comm_barrier() override;
  int replica_comm_recv(char* msg_data, int buf_len, int src_rep) override;
  int replica_comm_send(char* msg_data, int msg_len, int dest_rep) override;

  int check_atom_name_selections_available() override;
  int init_atom(int atom_number) override;
  int check_atom_id(int atom_number) override;
  int init_atom(cvm::residue_id const &residue,
                std::string const     &atom_name,
                std::string const     &segment_id) override;
  int check_atom_id(cvm::residue_id const &residue,
                    std::string const     &atom_name,
                    std::string const     &segment_id) override;
  void clear_atom(int index) override;

  void update_atom_properties(int index);

  cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                 cvm::atom_pos const &pos2) const;

  int load_atoms_pdb(char const *filename,
                     cvm::atom_group &atoms,
                     std::string const &pdb_field,
                     double pdb_field_value) override;

  int load_coords_pdb(char const *filename,
                      std::vector<cvm::atom_pos> &pos,
                      const std::vector<int> &indices,
                      std::string const &pdb_field,
                      double const pdb_field_value) override;


  int scalable_group_coms() override
  {
    return COLVARS_OK;
  }
  int init_atom_group(std::vector<int> const &atoms_ids) override;
  void clear_atom_group(int index) override;

  int update_group_properties(int index);

#if NAMD_VERSION_NUMBER >= 34471681

  int check_volmaps_available() override;

  int request_engine_volmap_by_id(int volmap_id) override;

  int request_engine_volmap_by_name(std::string const &volmap_name) override;

  /// Add map to GlobalMaster client (if not already in them)
  void request_globalmaster_volmap(int volmap_id);

  int init_internal_volmap_by_id(int volmap_id) override;

  int init_internal_volmap_by_name(std::string const &volmap_name) override;

  int load_internal_volmap_from_file(std::string const &filename) override;

  void clear_volmap(int index) override;

  int compute_volmap(int flags,
                     int index,
                     cvm::atom_group* ag,
                     cvm::real *value,
                     cvm::real *atom_field) override;

  /// Abstraction of the two types of NAMD volumetric maps
  template<class T>
  void getGridForceGridValue(int flags,
                             T const *grid,
                             cvm::atom_group* ag,
                             cvm::real *value,
                             cvm::real *atom_field);

  /// Implementation of inner loop; allows for atom list computation and use
  template<class T, int flags>
  void GridForceGridLoop(T const *g,
                         cvm::atom_group* ag,
                         cvm::real *value,
                         cvm::real *atom_field);

#endif

  std::ostream &output_stream(std::string const &output_name,
                              std::string const description) override;

  int flush_output_stream(std::string const &output_name) override;

  int flush_output_streams() override;

  int close_output_stream(std::string const &output_name) override;

  int close_output_streams() override;

  int backup_file(char const *filename) override;

  /// Get value of alchemical lambda parameter from back-end
  int get_alch_lambda(cvm::real* lambda);

  /// Set value of alchemical lambda parameter in back-end
  int send_alch_lambda(void);

  /// Request energy computation every freq steps
  int request_alch_energy_freq(int const freq);

  /// Get energy derivative with respect to lambda
  int get_dE_dlambda(cvm::real* dE_dlambda);

protected:

  /// Pointers to internally managed maps (set to nullptr for maps loaded by NAMD)
  std::vector<std::unique_ptr<GridforceFullMainGrid>> internal_gridforce_grids_;
};


#endif
