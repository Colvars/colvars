// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_NAMD_H
#define COLVARPROXY_NAMD_H

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
  SimParameters const *simparams;

  /// Self-explained
  BigReal thermostat_temperature;

  /// NAMD-style PRNG object
  Random random;

  /// How often NAMD is instructed to write state files
  size_t restart_frequency_s;

  bool first_timestep;
  size_t previous_NAMD_step;

  /// Used to submit restraint energy as MISC
  SubmitReduction *reduction;

  void init_tcl_pointers();

public:

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
  virtual int replica_num();
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

  std::ostream * output_stream(std::string const &output_name,
                               std::ios_base::openmode mode);
  int flush_output_stream(std::ostream *os);
  int close_output_stream(std::string const &output_name);
  int backup_file(char const *filename);

  char const *script_obj_to_str(unsigned char *obj);
  std::vector<std::string> script_obj_to_str_vector(unsigned char *obj);
};


#endif
