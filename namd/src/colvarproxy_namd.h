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

#include "colvarproxy_namd_version.h"

#include "Vector.h"
#include "ResizeArray.h"
#include "NamdTypes.h"
#include "SimParameters.h"
#include "Lattice.h"
#include "GlobalMaster.h"
#include "Random.h"
#include "ConfigList.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"

#define GLOBAL_MASTER_CKLOOP_CALC_ITEM 2000
#define GLOBAL_MASTER_CKLOOP_CALC_BIASES 2001
#define GLOBAL_MASTER_CKLOOP_CALC_SCRIPTED_BIASES 2002

/// \brief Communication between colvars and NAMD (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_namd : public colvarproxy, public GlobalMaster {

protected:

  /// \brief Array of atom indices (relative to the colvarproxy arrays),
  /// usedfor faster copy of atomic data
  std::vector<int> atoms_map;

  /// Pointer to the NAMD simulation input object
  SimParameters *simparams;

  /// Pointer to Controller object
  Controller const *controller;

  /// NAMD-style PRNG object
  Random random;

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

  friend class cvm::atom;

  colvarproxy_namd();
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

  cvm::real rand_gaussian() override
  {
    return random.gaussian();
  }

  cvm::real get_accelMD_factor() const override {
    return amd_weight_factor;
  }

  bool accelMD_enabled() const override {
    return accelMDOn;
  }

#if CMK_SMP && USE_CKLOOP
  colvarproxy::smp_mode_t get_smp_mode() const override;

  int set_smp_mode(smp_mode_t mode) override;

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
#ifdef COLVARS_USE_SOA
  int load_atoms_pdb(char const *filename,
                     cvm::atom_group_soa &atoms,
                     std::string const &pdb_field,
                     double pdb_field_value) override;
#else
  int load_atoms_pdb(char const *filename,
                     cvm::atom_group &atoms,
                     std::string const &pdb_field,
                     double const pdb_field_value) override;
#endif // COLVARS_USE_SOA


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

  int init_volmap_by_id(int volmap_id) override;

  int init_volmap_by_name(const char *volmap_name) override;

  int check_volmap_by_id(int volmap_id) override;

  int check_volmap_by_name(char const *volmap_name) override;

  int get_volmap_id_from_name(char const *volmap_name) override;

  void clear_volmap(int index) override;

#ifdef COLVARS_USE_SOA
  int compute_volmap(int flags,
                     int volmap_id,
                     cvm::atom_group_soa* ag,
                     cvm::real *value,
                     cvm::real *atom_field) override;
#else
  int compute_volmap(int flags,
                     int volmap_id,
                     cvm::atom_iter atom_begin,
                     cvm::atom_iter atom_end,
                     cvm::real *value,
                     cvm::real *atom_field) override;
#endif // COLVARS_USE_SOA

  /// Abstraction of the two types of NAMD volumetric maps
  template<class T>
  void getGridForceGridValue(int flags,
                             T const *grid,
#ifdef COLVARS_USE_SOA
                             cvm::atom_group_soa* ag,
#else
                             cvm::atom_iter atom_begin,
                             cvm::atom_iter atom_end,
#endif // COLVARS_USE_SOA
                             cvm::real *value,
                             cvm::real *atom_field);

  /// Implementation of inner loop; allows for atom list computation and use
  template<class T, int flags>
  void GridForceGridLoop(T const *g,
#ifdef COLVARS_USE_SOA
                         cvm::atom_group_soa* ag,
#else
                         cvm::atom_iter atom_begin,
                         cvm::atom_iter atom_end,
#endif // COLVARS_USE_SOA
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

};


#endif
