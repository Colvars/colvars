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

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"

/// \brief Communication between colvars and NAMD (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_namd : public colvarproxy, public GlobalMaster {

protected:

  /// \brief Array of atom indices (relative to the colvarproxy arrays),
  /// usedfor faster copy of atomic data
  std::vector<int> atoms_map;

  /// Pointer to the NAMD simulation input object
  SimParameters *simparams;

  /// Self-explained
  BigReal thermostat_temperature;

  /// NAMD-style PRNG object
  Random random;

  bool first_timestep;
  size_t previous_NAMD_step;

  /// Used to submit restraint energy as MISC
  SubmitReduction *reduction;

  /// Accelerated MD reweighting factor
  bool accelMDOn;
  cvm::real amd_weight_factor;
  void update_accelMD_info();

public:

  virtual void init_tcl_pointers();

  friend class cvm::atom;

  colvarproxy_namd();
  ~colvarproxy_namd();

  int setup();
  int reset();

  // synchronize the local arrays with requested or forced atoms
  int update_atoms_map(AtomIDList::const_iterator begin,
                       AtomIDList::const_iterator end);

  void calculate();

  void log(std::string const &message);
  void error(std::string const &message);
  int set_unit_system(std::string const &units_in, bool check_only);
  void exit(std::string const &message);
  void add_energy(cvm::real energy);
  void request_total_force(bool yesno);

  bool total_forces_enabled() const
  {
    return total_force_requested;
  }

  int run_force_callback();
  int run_colvar_callback(std::string const &name,
                          std::vector<const colvarvalue *> const &cvcs,
                          colvarvalue &value);
  int run_colvar_gradient_callback(std::string const &name,
                                   std::vector<const colvarvalue *> const &cvcs,
                                   std::vector<cvm::matrix2d<cvm::real> > &gradient);

  cvm::real backend_angstrom_value()
  {
    return 1.0;
  }

  cvm::real boltzmann()
  {
    return 0.001987191;
  }

  cvm::real temperature()
  {
    return thermostat_temperature;
  }

  cvm::real rand_gaussian()
  {
    return random.gaussian();
  }

  cvm::real dt()
  {
    return simparams->dt;
  }

  virtual cvm::real get_accelMD_factor() const {
    return amd_weight_factor;
  }

  bool accelMD_enabled() const {
    return accelMDOn;
  }

#if CMK_SMP && USE_CKLOOP
  int smp_enabled()
  {
    if (b_smp_active) {
      return COLVARS_OK;
    }
    return COLVARS_ERROR;
  }

  int smp_colvars_loop();

  int smp_biases_loop();

  int smp_biases_script_loop();

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
    charm_lock_state = CmiCreateLock();
    return COLVARS_OK;
  }

  int smp_trylock()
  {
    return COLVARS_NOT_IMPLEMENTED;
  }

  int smp_unlock()
  {
    CmiDestroyLock(charm_lock_state);
    return COLVARS_OK;
  }

#endif // #if CMK_SMP && USE_CKLOOP

  virtual int replica_enabled();
  virtual int replica_index();
  virtual int num_replicas();
  virtual void replica_comm_barrier();
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep);
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep);

  int init_atom(int atom_number);
  int check_atom_id(int atom_number);
  int init_atom(cvm::residue_id const &residue,
                std::string const     &atom_name,
                std::string const     &segment_id);
  int check_atom_id(cvm::residue_id const &residue,
                    std::string const     &atom_name,
                    std::string const     &segment_id);
  void clear_atom(int index);

  void update_atom_properties(int index);

  cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                 cvm::atom_pos const &pos2) const;

  int load_atoms(char const *filename,
                 cvm::atom_group &atoms,
                 std::string const &pdb_field,
                 double const pdb_field_value = 0.0);

  int load_coords(char const *filename,
                  std::vector<cvm::atom_pos> &pos,
                  const std::vector<int> &indices,
                  std::string const &pdb_field,
                  double const pdb_field_value = 0.0);


  int scalable_group_coms()
  {
    return COLVARS_OK;
  }
  int init_atom_group(std::vector<int> const &atoms_ids);
  void clear_atom_group(int index);
  int update_group_properties(int index);

#if NAMD_VERSION_NUMBER >= 34471681

  virtual int init_volmap_by_id(int volmap_id);

  virtual int init_volmap_by_name(const char *volmap_name);

  virtual int check_volmap_by_id(int volmap_id);

  virtual int check_volmap_by_name(char const *volmap_name);

  virtual int get_volmap_id_from_name(char const *volmap_name);

  virtual void clear_volmap(int index);

  virtual int compute_volmap(int flags,
                             int volmap_id,
                             cvm::atom_iter atom_begin,
                             cvm::atom_iter atom_end,
                             cvm::real *value,
                             cvm::real *atom_field);

  /// Abstraction of the two types of NAMD volumetric maps
  template<class T>
  void getGridForceGridValue(int flags,
                             T const *grid,
                             cvm::atom_iter atom_begin,
                             cvm::atom_iter atom_end,
                             cvm::real *value,
                             cvm::real *atom_field);

  /// Implementation of inner loop; allows for atom list computation and use
  template<class T, int flags>
  void GridForceGridLoop(T const *g,
                         cvm::atom_iter atom_begin,
                         cvm::atom_iter atom_end,
                         cvm::real *value,
                         cvm::real *atom_field);

#endif

  std::ostream * output_stream(std::string const &output_name,
                               std::ios_base::openmode mode);
  int flush_output_stream(std::ostream *os);
  int close_output_stream(std::string const &output_name);
  int backup_file(char const *filename);

};


#endif
