// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_H
#define COLVARPROXY_H

#include <functional>
#include <list>
#include <map>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvar_gpu_support.h"

/// \file colvarproxy.h
/// \brief Colvars proxy classes
///
/// This file declares the class for the object responsible for interfacing
/// Colvars with other codes (MD engines, VMD, Python).  The \link colvarproxy
/// \endlink class is a derivative of multiple classes, each devoted to a
/// specific task (e.g. \link colvarproxy_atoms \endlink to access data for
/// individual atoms).
///
/// To interface to a new MD engine, the simplest solution is to derive a new
/// class from \link colvarproxy \endlink.  Currently implemented are: \link
/// colvarproxy_lammps, \endlink, \link colvarproxy_namd, \endlink, \link
/// colvarproxy_vmd \endlink.


// forward declarations
class colvarscript;

#if defined(_OPENMP)
#include <omp.h>
#else
struct omp_lock_t;
#endif

#ifdef COLVARS_MPI
#include <mpi.h>
typedef MPI_Comm replicas_mpi_comm_t;
#else
typedef void * replicas_mpi_comm_t;
#endif

#if defined(NAMD_TCL) || defined(VMDTCL)
#define COLVARS_TCL
#endif

#ifdef COLVARS_TCL
#include <tcl.h>
#else
// Allow for placeholders Tcl_Interp* variables
typedef void Tcl_Interp;
#endif


/// Interface between Colvars and MD engine (GROMACS, LAMMPS, NAMD, VMD...)
///
/// This is the base class: each engine is supported by a derived class.
class colvarproxy
{

public:

  /// \brief Name of the unit system used internally by Colvars (by default, that of the back-end).
  /// Supported depending on the back-end: real (A, kcal/mol), metal (A, eV), electron (Bohr, Hartree), gromacs (nm, kJ/mol)
  /// Note: calls to back-end PBC functions assume back-end length unit
  /// We use different unit from back-end in VMD bc using PBC functions from colvarproxy base class
  /// Colvars internal units are user specified, because the module exchanges info in unknown
  /// composite dimensions with user input, while it only exchanges quantities of known
  /// dimension with the back-end (length and forces)
  std::string units;

  /// \brief Request to set the units used internally by Colvars
  virtual int set_unit_system(std::string const &units, bool check_only);

  /// \brief Convert a length from Angstrom to internal
  inline cvm::real angstrom_to_internal(cvm::real l) const
  {
    return l * angstrom_value_;
  }

  /// \brief Convert a length from internal to Angstrom
  inline cvm::real internal_to_angstrom(cvm::real l) const
  {
    return l / angstrom_value_;
  }

  /// Boltzmann constant, with unit the same as energy / K
  inline cvm::real boltzmann() const
  {
    return boltzmann_;
  }

  /// Current target temperature of the simulation (K units)
  inline cvm::real target_temperature() const
  {
    return target_temperature_;
  }

  /// Set the current target temperature of the simulation (K units)
  virtual int set_target_temperature(cvm::real T);

  /// Time step of the simulation (fs units)
  inline double dt() const
  {
    return timestep_;
  }

  /// Set the current integration timestep of the simulation (fs units)
  virtual int set_integration_timestep(cvm::real dt);

  /// Time step of the simulation (fs units)
  inline int time_step_factor() const
  {
    return time_step_factor_;
  }

  /// Set the current integration timestep of the simulation (fs units)
  virtual int set_time_step_factor(int fact);

  /// \brief Pseudo-random number with Gaussian distribution
  virtual cvm::real rand_gaussian(void);

  /// Pass restraint energy value for current timestep to MD engine
  virtual void add_energy(cvm::real energy);

  /// Use the PBC functions from the Colvars library (as opposed to MD engine)
  inline bool & use_internal_pbc() { return use_internal_pbc_; }

  /// \brief Get the PBC-aware distance vector between two positions
  virtual cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                         cvm::atom_pos const &pos2) const;

  /// Inline version of position_distance()
  cvm::rvector position_distance_internal(cvm::atom_pos const &pos1,
                                          cvm::atom_pos const &pos2) const;

  /// Kernel used by position_distance_internal()
  /// @param pos1 First position
  /// @param pos2 Second position
  /// @param a Unit cell vector
  /// @param b Unit cell vector
  /// @param c Unit cell vector
  /// @param a_r Reciprocal cell vector
  /// @param b_r Reciprocal cell vector
  /// @param c_r Reciprocal cell vector
  /// @param a_p Is this dimension periodic?
  /// @param b_p Is this dimension periodic?
  /// @param c_p Is this dimension periodic?
  /// @return (pos2 - pos1) without periodicity, minimum-image difference otherwise
  static cvm::rvector position_distance_kernel(cvm::atom_pos const &pos1,
                                               cvm::atom_pos const &pos2,
                                               cvm::rvector const &a,
                                               cvm::rvector const &b,
                                               cvm::rvector const &c,
                                               cvm::rvector const &a_r,
                                               cvm::rvector const &b_r,
                                               cvm::rvector const &c_r,
                                               bool a_p = false,
                                               bool b_p = false,
                                               bool c_p = false);


  /// Recompute PBC reciprocal lattice (assumes XYZ periodicity)
  void update_pbc_lattice();

  /// Set the lattice vectors to zero
  void reset_pbc_lattice();

  /// \brief Tell the proxy whether total forces are needed (they may not
  /// always be available)
  virtual void request_total_force(bool yesno);

  /// Are total forces being used?
  virtual bool total_forces_enabled() const;

  /// Are total forces from the current step available?
  /// in which case they are really system forces
  virtual bool total_forces_same_step() const;

  /// Get the molecule ID when called in VMD; raise error otherwise
  /// \param molid Set this argument equal to the current VMD molid
  virtual int get_molid(int &molid);

  /// Get value of alchemical lambda parameter from back-end (if available)
  virtual int get_alch_lambda(cvm::real* lambda);

  /// Set value of alchemical lambda parameter to be sent to back-end at end of timestep
  void set_alch_lambda(cvm::real lambda);

  /// Send cached value of alchemical lambda parameter to back-end (if available)
  virtual int send_alch_lambda();

  /// Request energy computation every freq steps (necessary for NAMD3, not all back-ends)
  virtual int request_alch_energy_freq(int const /* freq */) {
    return COLVARS_OK;
  }

  /// Get energy derivative with respect to lambda (if available)
  virtual int get_dE_dlambda(cvm::real* dE_dlambda);

  /// Apply a scalar force on dE_dlambda (back-end distributes it onto atoms)
  virtual int apply_force_dE_dlambda(cvm::real* force);

  /// Get energy second derivative with respect to lambda (if available)
  virtual int get_d2E_dlambda2(cvm::real* d2E_dlambda2);

  /// Force to be applied onto alch. lambda, propagated from biasing forces on dE_dlambda
  cvm::real indirect_lambda_biasing_force;

  /// Get weight factor from accelMD
  virtual cvm::real get_accelMD_factor() const {
    cvm::error_static(cvmodule, "Error: accessing the reweighting factor of accelerated MD  "
               "is not yet implemented in the MD engine.\n",
               COLVARS_NOT_IMPLEMENTED);
    return 1.0;
  }
  virtual bool accelMD_enabled() const {
    return false;
  }

#if ( defined(COLVARS_CUDA) || defined(COLVARS_HIP) )
  template <typename T>
  using allocator_type = colvars_gpu::CudaHostAllocator<T>;
#else
  template <typename T>
  using allocator_type = std::allocator<T>;
#endif
  using atom_buffer_real_t = std::vector<cvm::real, allocator_type<cvm::real>>;
  using atom_buffer_rvector_t = std::vector<cvm::rvector, allocator_type<cvm::rvector>>;

  /// Prepare this atom for collective variables calculation, selecting it by
  /// numeric index (1-based)
  virtual int init_atom(int atom_number);

  /// Check that this atom number is valid, but do not initialize the
  /// corresponding atom yet
  virtual int check_atom_id(int atom_number);

  /// Check whether it is possible to select atoms by residue number name
  virtual int check_atom_name_selections_available();

  /// Select this atom for collective variables calculation, using name and
  /// residue number.  Not all programs support this: leave this function as
  /// is in those cases.
  virtual int init_atom(cvm::residue_id const &residue,
                        std::string const     &atom_name,
                        std::string const     &segment_id);

  /// Check that this atom is valid, but do not initialize it yet
  virtual int check_atom_id(cvm::residue_id const &residue,
                            std::string const     &atom_name,
                            std::string const     &segment_id);

  /// \brief Used by the atom class destructor: rather than deleting the array slot
  /// (costly) set the corresponding atoms_refcount to zero
  virtual void clear_atom(int index);

  /// Clear atomic data
  int reset_atoms();

  /// Get the numeric ID of the given atom
  /// \param index Internal index in the Colvars arrays
  inline int get_atom_id(int index) const
  {
    return atoms_ids[index];
  }

  /// Get the mass of the given atom
  /// \param index Internal index in the Colvars arrays
  inline cvm::real get_atom_mass(int index) const
  {
    return atoms_masses[index];
  }

  /// Increase the reference count of the given atom
  /// \param index Internal index in the Colvars arrays
  inline void increase_refcount(int index)
  {
    atoms_refcount[index] += 1;
  }

  /// Get the charge of the given atom
  /// \param index Internal index in the Colvars arrays
  inline cvm::real get_atom_charge(int index) const
  {
    return atoms_charges[index];
  }

  /// Read the current position of the given atom
  /// \param index Internal index in the Colvars arrays
  inline cvm::rvector get_atom_position(int index) const
  {
    return atoms_positions[index];
  }

  /// Read the current total force of the given atom
  /// \param index Internal index in the Colvars arrays
  inline cvm::rvector get_atom_total_force(int index) const
  {
    return atoms_total_forces[index];
  }

  /// Request that this force is applied to the given atom
  /// \param index Internal index in the Colvars arrays
  /// \param new_force Force to add
  inline void apply_atom_force(int index, cvm::rvector const &new_force)
  {
    atoms_new_colvar_forces[index] += new_force;
  }

  /// Read the current velocity of the given atom
  inline cvm::rvector get_atom_velocity(int /* index */)
  {
    cvm::error_static(cvmodule, "Error: reading the current velocity of an atom "
               "is not yet implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
    return cvm::rvector(0.0);
  }

  inline std::vector<int> const *get_atom_ids() const
  {
    return &atoms_ids;
  }

  /// Return number of atoms with positive reference count
  size_t get_num_active_atoms() const;

  inline atom_buffer_real_t const *get_atom_masses() const
  {
    return &atoms_masses;
  }

  inline atom_buffer_real_t *modify_atom_masses()
  {
    // assume that we are requesting masses to change them
    updated_masses_ = true;
    return &atoms_masses;
  }

  inline atom_buffer_real_t const *get_atom_charges()
  {
    return &atoms_charges;
  }

  inline atom_buffer_real_t *modify_atom_charges()
  {
    // assume that we are requesting charges to change them
    updated_charges_ = true;
    return &atoms_charges;
  }

  inline atom_buffer_rvector_t const *get_atom_positions() const
  {
    return &atoms_positions;
  }

  inline atom_buffer_rvector_t *modify_atom_positions()
  {
    return &atoms_positions;
  }

  inline atom_buffer_rvector_t const *get_atom_total_forces() const
  {
    return &atoms_total_forces;
  }

  inline atom_buffer_rvector_t *modify_atom_total_forces()
  {
    return &atoms_total_forces;
  }

  inline atom_buffer_rvector_t const *get_atom_applied_forces() const
  {
    return &atoms_new_colvar_forces;
  }

  inline atom_buffer_rvector_t *modify_atom_applied_forces()
  {
    return &atoms_new_colvar_forces;
  }

  /// Compute the root-mean-square of the applied forces
  void compute_rms_atoms_applied_force();

  /// Compute the maximum norm among all applied forces
  void compute_max_atoms_applied_force();

  /// Get the root-mean-square of the applied forces
  inline cvm::real rms_atoms_applied_force() const
  {
    return atoms_rms_applied_force_;
  }

  /// Get the maximum norm among all applied forces
  inline cvm::real max_atoms_applied_force() const
  {
    return atoms_max_applied_force_;
  }

  /// Get the atom ID with the largest applied force
  inline int max_atoms_applied_force_id() const
  {
    return atoms_max_applied_force_id_;
  }

  /// Whether the atom list has been modified internally
  inline bool modified_atom_list() const
  {
    return modified_atom_list_;
  }

  /// Reset the modified atom list flag
  inline void reset_modified_atom_list()
  {
    modified_atom_list_ = false;
  }

  /// Record whether masses have been updated
  inline bool updated_masses() const
  {
    return updated_masses_;
  }

  /// Record whether masses have been updated
  inline bool updated_charges() const
  {
    return updated_charges_;
  }

  /// Clear atom group data
  int reset_atom_groups();

  /// \brief Whether this proxy implementation has capability for scalable groups
  virtual int scalable_group_coms();

  /// Prepare this group for collective variables calculation, selecting atoms by internal ids (0-based)
  virtual int init_atom_group(std::vector<int> const &atoms_ids);

  /// \brief Used by the atom_group class destructor
  virtual void clear_atom_group(int index);

  /// Get the numeric ID of the given atom group (for the MD program)
  inline int get_atom_group_id(int index) const
  {
    return atom_groups_ids[index];
  }

  /// Get the mass of the given atom group
  inline cvm::real get_atom_group_mass(int index) const
  {
    return atom_groups_masses[index];
  }

  /// Get the charge of the given atom group
  inline cvm::real get_atom_group_charge(int index) const
  {
    return atom_groups_charges[index];
  }

  /// Read the current position of the center of mass given atom group
  inline cvm::rvector get_atom_group_com(int index) const
  {
    return atom_groups_coms[index];
  }

  /// Read the current total force of the given atom group
  inline cvm::rvector get_atom_group_total_force(int index) const
  {
    return atom_groups_total_forces[index];
  }

  /// Request that this force is applied to the given atom group
  inline void apply_atom_group_force(int index, cvm::rvector const &new_force)
  {
    atom_groups_new_colvar_forces[index] += new_force;
  }

  /// Read the current velocity of the given atom group
  inline cvm::rvector get_atom_group_velocity(int /* index */)
  {
    cvm::error_static(cvmodule, "Error: reading the current velocity of an atom group is not yet implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
    return cvm::rvector(0.0);
  }

  inline std::vector<int> const *get_atom_group_ids() const
  {
    return &atom_groups_ids;
  }

  /// Return number of atom groups with positive reference count
  size_t get_num_active_atom_groups() const;

  inline std::vector<cvm::real> *modify_atom_group_masses()
  {
    // TODO updated_masses
    return &atom_groups_masses;
  }

  inline std::vector<cvm::real> *modify_atom_group_charges()
  {
    // TODO updated masses
    return &atom_groups_charges;
  }

  inline std::vector<cvm::rvector> *modify_atom_group_positions()
  {
    return &atom_groups_coms;
  }

  inline std::vector<cvm::rvector> *modify_atom_group_total_forces()
  {
    return &atom_groups_total_forces;
  }

  inline std::vector<cvm::rvector> *modify_atom_group_applied_forces()
  {
    return &atom_groups_new_colvar_forces;
  }

  /// Compute the root-mean-square of the applied forces
  void compute_rms_atom_groups_applied_force();

  /// Compute the maximum norm among all applied forces
  void compute_max_atom_groups_applied_force();

  /// Get the root-mean-square of the applied forces
  inline cvm::real rms_atom_groups_applied_force() const
  {
    return atom_groups_rms_applied_force_;
  }

  /// Get the maximum norm among all applied forces
  inline cvm::real max_atom_groups_applied_force() const
  {
    return atom_groups_max_applied_force_;
  }

  /// Clear volumetric map data
  int reset_volmaps();

  /// Test whether this implementation can use volumetric maps as CVs
  virtual int check_volmaps_available();

  /// Create a slot for a volumetric map not requested yet
  int add_volmap_slot(int volmap_id);

  /// Request and prepare this volumetric map for use by Colvars
  /// \param volmap_id Numeric ID used by the MD engine
  /// \returns Index of the map in the colvarproxy arrays
  virtual int init_volmap_by_id(int volmap_id);

  /// Request and prepare this volumetric map for use by Colvars
  /// \param volmap_name Name used by the MD engine
  /// \returns Index of the map in the colvarproxy arrays
  virtual int init_volmap_by_name(char const *volmap_name);

  /// Check that the given volmap ID is valid (return COLVARS_OK if it is)
  /// \param volmap_id Numeric ID used by the MD engine
  /// \returns Error code
  virtual int check_volmap_by_id(int volmap_id);

  /// Check that the given volmap name is valid (return COLVARS_OK if it is)
  /// \param volmap_name Name used by the MD engine
  /// \returns Error code
  virtual int check_volmap_by_name(char const *volmap_name);

  /// Request and prepare this volumetric map for use by Colvars
  int init_volmap_by_name(std::string const &volmap_name);

  /// Check that the given volmap name is valid (return COLVARS_OK if it is)
  int check_volmap_by_name(std::string const &volmap_name);

  /// \brief Used by the CVC destructors
  virtual void clear_volmap(int index);

  /// Get the numeric ID of the given volumetric map (for the MD program)
  virtual int get_volmap_id_from_name(char const *volmap_name);

  /// Get the numeric ID of the given volumetric map (for the MD program)
  inline int get_volmap_id(int index) const
  {
    return volmaps_ids[index];
  }

  /// Read the current value of the volumetric map
  inline cvm::real get_volmap_value(int index) const
  {
    return volmaps_values[index];
  }

  /// Request that this force is applied to the given volumetric map
  inline void apply_volmap_force(int index, cvm::real const &new_force)
  {
    volmaps_new_colvar_forces[index] += new_force;
  }

  /// Re-weigh an atomic field (e.g. a colvar) by the value of a volumetric map

  /// \param flags Combination of flags
  /// \param volmap_id Numeric index of the map (no need to request it)
  /// \param ag Pointer to the SOA atom group
  /// \param value Pointer to location of total to increment
  /// \param atom_field Array of atomic field values (if NULL, ones are used)
  virtual int compute_volmap(int flags,
                             int volmap_id,
                             cvm::atom_group* ag,
                             cvm::real *value,
                             cvm::real *atom_field);

  /// Flags controlling what computation is done on the map
  enum {
    volmap_flag_null = 0,
    volmap_flag_gradients = 1,
    volmap_flag_use_atom_field = (1<<8)
  };

  /// Compute the root-mean-square of the applied forces
  void compute_rms_volmaps_applied_force();

  /// Compute the maximum norm among all applied forces
  void compute_max_volmaps_applied_force();

  enum class smp_mode_t {cvcs, inner_loop, gpu, none};

  /// Get the current SMP mode
  virtual smp_mode_t get_smp_mode() const;

  /// Get available SMP modes
  virtual std::vector<smp_mode_t> get_available_smp_modes() const;

  /// Get the preferred SMP mode
  virtual smp_mode_t get_preferred_smp_mode() const;

  /// Set the current SMP mode
  virtual int set_smp_mode(smp_mode_t mode);

  /// Distribute computation over threads using OpenMP, unless overridden in the backend (e.g. NAMD)
  virtual int smp_loop(int n_items, std::function<int (int)> const &worker);

  /// Distribute calculation of biases across threads
  virtual int smp_biases_loop();

  /// Distribute calculation of biases across threads 2nd through last, with all scripted biased on 1st thread
  virtual int smp_biases_script_loop();

  /// Index of this thread
  virtual int smp_thread_id();

  /// Number of threads sharing this address space
  virtual int smp_num_threads();

  /// Lock the proxy's shared data for access by a thread, if threads are implemented; if not implemented, does nothing
  virtual int smp_lock();

  /// Attempt to lock the proxy's shared data
  virtual int smp_trylock();

  /// Release the lock
  virtual int smp_unlock();

  /// Set the multiple replicas communicator
  virtual void set_replicas_mpi_communicator(replicas_mpi_comm_t comm);

  /// Indicate if multi-replica support is available and active
  virtual int check_replicas_enabled();

  /// Index of this replica
  virtual int replica_index();

  /// Total number of replicas
  virtual int num_replicas();

  /// Synchronize replica with others
  virtual void replica_comm_barrier();

  /// Receive data from other replica
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep);

  /// Send data to other replica
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep);

  /// Pointer to the scripting interface object
  /// (does not need to be allocated in a new interface)
  colvarscript *script = nullptr;

  /// Run a user-defined colvar forces script
  virtual int run_force_callback();

  virtual int run_colvar_callback(
                std::string const &name,
                std::vector<const colvarvalue *> const &cvcs,
                colvarvalue &value);

  virtual int run_colvar_gradient_callback(
                std::string const &name,
                std::vector<const colvarvalue *> const &cvcs,
                std::vector<cvm::matrix2d<cvm::real> > &gradient);

  /// Is Tcl available? (trigger initialization if needed)
  inline bool tcl_available() {
#if defined(COLVARS_TCL)
    return true;
#else
    return false;
#endif
  }

  /// Get a string representation of the Tcl object pointed to by obj
  char const *tcl_get_str(void *obj);

  int tcl_run_script(std::string const &script);

  int tcl_run_file(std::string const &fileName);

  /// Tcl implementation of run_force_callback()
  int tcl_run_force_callback();

  /// Tcl implementation of run_colvar_callback()
  int tcl_run_colvar_callback(
              std::string const &name,
              std::vector<const colvarvalue *> const &cvcs,
              colvarvalue &value);

  /// Tcl implementation of run_colvar_gradient_callback()
  int tcl_run_colvar_gradient_callback(
              std::string const &name,
              std::vector<const colvarvalue *> const &cvcs,
              std::vector<cvm::matrix2d<cvm::real> > &gradient);

  /// Get a pointer to the Tcl interpreter
  inline Tcl_Interp *get_tcl_interp()
  {
    return tcl_interp_;
  }

  /// Set the pointer to the Tcl interpreter
  inline void set_tcl_interp(Tcl_Interp *interp)
  {
    tcl_interp_ = interp;
  }

  /// Set Tcl pointers
  virtual void init_tcl_pointers();

  /// Ensure that we're on the main thread (derived class will do actual check)
  virtual bool io_available();

  /// \brief Save the current frame number in the argument given
  // Returns error code
  virtual int get_frame(long int &);

  /// \brief Set the current frame number (as well as cvmodule->it)
  // Returns error code
  virtual int set_frame(long int);

  /// Get the current working directory of this process
  std::string get_current_work_dir() const;

  /// Join two paths using the operating system's path separation
  std::string join_paths(std::string const &path1, std::string const &path2) const;

  /// \brief Rename the given file, before overwriting it
  virtual int backup_file(char const *filename);

  /// \brief Rename the given file, before overwriting it
  inline int backup_file(std::string const &filename)
  {
    return backup_file(filename.c_str());
  }

  /// Remove the given file (on Windows only, rename to filename.old)
  virtual int remove_file(char const *filename);

  /// Remove the given file (on Windows only, rename to filename.old)
  inline int remove_file(std::string const &filename)
  {
    return remove_file(filename.c_str());
  }

  /// Rename the given file
  virtual int rename_file(char const *filename, char const *newfilename);

  /// Rename the given file
  inline int rename_file(std::string const &filename,
                         std::string const &newfilename)
  {
    return rename_file(filename.c_str(), newfilename.c_str());
  }

  /// Prefix of the input state file to be read next
  inline std::string const & input_prefix() const
  {
    return input_prefix_str;
  }

  /// Initialize input_prefix (NOTE: it will be erased after state file is read)
  virtual int set_input_prefix(std::string const &prefix);

  /// Default prefix to be used for all output files (final configuration)
  inline std::string const & output_prefix() const
  {
    return output_prefix_str;
  }

  /// Set default output prefix
  virtual int set_output_prefix(std::string const &prefix);

  /// Prefix of the restart (checkpoint) file to be written next
  inline std::string const & restart_output_prefix() const
  {
    return restart_output_prefix_str;
  }

  /// Set default restart state file prefix
  virtual int set_restart_output_prefix(std::string const &prefix);

  /// Default restart frequency (as set by the simulation engine)
  inline int default_restart_frequency() const
  {
    return restart_frequency_engine;
  }

  /// Communicate/set the restart frequency of the simulation engine
  virtual int set_default_restart_frequency(int freq);

  // The input stream functions below are not virtual, because they currently
  // rely on the fact that either std::ifstream or std::istringstream is used

  /// Returns a reference to given input stream, creating it if needed
  /// \param input_name File name (later only a handle)
  /// \param description Purpose of the file
  /// \param error_on_fail Raise error when failing to open (allow testing)
  std::istream & input_stream(std::string const &input_name,
                              std::string const description = "file/channel",
                              bool error_on_fail = true);

  /// Returns a reference to given input stream, creating it if needed
  /// \param input_name Identifier of the input stream
  /// \param content Set this string as the stream's buffer
  /// \param description Purpose of the stream
  std::istream & input_stream_from_string(std::string const &input_name,
                                          std::string const &content,
                                          std::string const description = "string");

  /// Check if the file/channel is open (without opening it if not)
  bool input_stream_exists(std::string const &input_name);

  /// Closes the given input stream
  int close_input_stream(std::string const &input_name);

  /// Closes all input streams
  int close_input_streams();

  /// Same as close_input_stream(), but also removes the corresponding entry from memory
  int delete_input_stream(std::string const &input_name);

  /// List all input streams that were opened at some point
  std::list<std::string> list_input_stream_names() const;

  /// Returns a reference to the named output file/channel (open it if needed)
  /// \param output_name File name or identifier
  /// \param description Purpose of the file
  virtual std::ostream &output_stream(std::string const &output_name,
                                      std::string const description);

  /// Check if the file/channel is open (without opening it if not)
  virtual bool output_stream_exists(std::string const &output_name);

  /// Flushes the given output file/channel
  virtual int flush_output_stream(std::string const &output_name);

  /// Flushes all output files/channels
  virtual int flush_output_streams();

  /// Closes the given output file/channel
  virtual int close_output_stream(std::string const &output_name);

  /// Close all open files/channels to prevent data loss
  virtual int close_output_streams();

  /// \brief Whether the proxy supports GPU
  bool has_gpu_support() const {
    return support_gpu;
  }
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP) || defined (COLVARS_SYCL)
  /// \brief Get the default CUDA stream from the proxy
  virtual cudaStream_t get_default_stream() {return (cudaStream_t)0;}
  /**
   * @brief Template function to allocate host-pinned memory
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated host-pinned memory
   * @param[in] len Number of elements to allocate
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int allocate_host(T **pp, const size_t len) {
    return allocate_host_T((void **)pp, len, sizeof(T));
  }
  /**
   * @brief Template function to deallocate host-pinned memory
   *
   * @tparam T The type of elements to deallocate
   * @param[in,out] pp Pointer to the pointer that holds the allocated host-pinned memory
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int deallocate_host(T **pp) {
    return deallocate_host_T((void **)pp);
  }
  /**
   * @brief Template function to allocate device memory
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated device memory
   * @param[in] len Number of elements to allocate
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int allocate_device(T **pp, const size_t len) {
    return allocate_device_T((void **)pp, len, sizeof(T));
  }
  /**
   * @brief Template function to reallocate device memory
   *
   * This function first deallocates any existing memory pointed to by `*pp`,
   * then allocates new device memory for `len` elements of type `T`.
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated device memory
   * @param[in] len Number of elements to allocate
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int reallocate_device(T **pp, const size_t len) {
    int error_code = COLVARS_OK;
    error_code |= deallocate_device(pp);
    error_code |= allocate_device_T((void **)pp, len, sizeof(T));
    return error_code;
  }
  /**
   * @brief Template function to reallocate host-pinned memory
   *
   * This function first deallocates any existing memory pointed to by `*pp`,
   * then allocates new host-pinned memory for `len` elements of type `T`.
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated host-pinned memory
   * @param[in] len Number of elements to allocate
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int reallocate_host(T **pp, const size_t len) {
    int error_code = COLVARS_OK;
    error_code |= deallocate_host(pp);
    error_code |= allocate_host_T((void **)pp, len, sizeof(T));
    return error_code;
  }
  /**
   * @brief Template function to allocate device memory asynchronously
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated device memory
   * @param[in] len Number of elements to allocate
   * @param[in] stream The CUDA stream to use for the allocation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int allocate_device_async(T **pp, const size_t len, cudaStream_t stream) {
    return allocate_device_T_async((void **)pp, len, sizeof(T), stream);
  }
  /**
   * @brief Template function to deallocate device memory
   *
   * @tparam T The type of elements to deallocate
   * @param[in,out] pp Pointer to the pointer that holds the allocated device memory
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int deallocate_device(T **pp) {
    return deallocate_device_T((void **)pp);
  }
  /**
   * @brief Template function to deallocate device memory asynchronously
   *
   * @tparam T The type of elements to deallocate
   * @param[in,out] pp Pointer to the pointer that holds the allocated device memory
   * @param[in] stream The CUDA stream to use for the deallocation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int deallocate_device_async(T **pp, cudaStream_t stream) {
    return deallocate_device_T_async((void **)pp, stream);
  }
  /**
   * @brief Template function to clear a device array to zero
   *
   * @tparam T The type of elements in the array
   * @param[in] data Pointer to the device array to clear
   * @param[in] ndata Number of elements in the array
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int clear_device_array(T *data, const size_t ndata) {
    return clear_device_array_T(data, ndata, sizeof(T));
  }
  /**
   * @brief Template function to clear a device array to zero asynchronously
   *
   * @tparam T The type of elements in the array
   * @param[in] data Pointer to the device array to clear
   * @param[in] ndata Number of elements in the array
   * @param[in] stream The CUDA stream to use for the operation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int clear_device_array_async(T *data, const size_t ndata, cudaStream_t stream) {
    return clear_device_array_T_async(data, ndata, sizeof(T), stream);
  }
  /**
   * @brief Template function to copy data from host to device
   *
   * @tparam T The type of elements to copy
   * @param[in] h_array Pointer to the host array
   * @param[in] d_array Pointer to the device array
   * @param[in] array_len Number of elements to copy
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_HtoD(const T *h_array, T *d_array, size_t array_len) {
    return copy_HtoD_T(h_array, d_array, array_len, sizeof(T));
  }
  /**
   * @brief Template function to copy data from host to device asynchronously
   *
   * @tparam T The type of elements to copy
   * @param[in] h_array Pointer to the host array
   * @param[out] d_array Pointer to the device array
   * @param[in] array_len Number of elements to copy
   * @param[in] stream The CUDA stream to use for the operation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_HtoD_async(const T *h_array, T *d_array, size_t array_len, cudaStream_t stream) {
    return copy_HtoD_T_async(h_array, d_array, array_len, sizeof(T), stream);
  }
  /**
   * @brief Template function to copy data from device to host
   *
   * @tparam T The type of elements to copy
   * @param[in] d_array Pointer to the device array
   * @param[out] h_array Pointer to the host array
   * @param[in] array_len Number of elements to copy
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_DtoH(const T *d_array, T *h_array, size_t array_len) {
    return copy_DtoH_T(d_array, h_array, array_len, sizeof(T));
  }
  /**
   * @brief Template function to copy data from device to host asynchronously
   *
   * @tparam T The type of elements to copy
   * @param[in] d_array Pointer to the device array
   * @param[out] h_array Pointer to the host array
   * @param[in] array_len Number of elements to copy
   * @param[in] stream The CUDA stream to use for the operation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_DtoH_async(const T *d_array, T *h_array, size_t array_len, cudaStream_t stream) {
    return copy_DtoH_T_async(d_array, h_array, array_len, sizeof(T), stream);
  }
  /**
   * @brief Template function to copy data from device to device
   *
   * @tparam T The type of elements to copy
   * @param[in] d_src Pointer to the source device array
   * @param[out] d_dst Pointer to the destination device array
   * @param[in] array_len Number of elements to copy
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_DtoD(const T *d_src, T *d_dst, size_t array_len) {
    return copy_DtoD_T(d_src, d_dst, array_len, sizeof(T));
  }
  /**
   * @brief Template function to copy data from device to device asynchronously
   *
   * @tparam T The type of elements to copy
   * @param[in] d_src Pointer to the source device array
   * @param[out] d_dst Pointer to the destination device array
   * @param[in] array_len Number of elements to copy
   * @param[in] stream The CUDA stream to use for the operation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_DtoD_async(const T *d_src, T *d_dst, size_t array_len, cudaStream_t stream) {
    return copy_DtoD_T_async(d_src, d_dst, array_len, sizeof(T), stream);
  }
  /// @brief Memory management and data transfer implementations
  /// @{
  virtual int allocate_host_T(void **pp, const size_t len, const size_t sizeofT);
  virtual int deallocate_host_T(void **pp);
  virtual int allocate_device_T(void **pp, const size_t len, const size_t sizeofT);
  virtual int deallocate_device_T(void **pp);
  virtual int clear_device_array_T(void *data, const size_t ndata, const size_t sizeofT);
  virtual int allocate_device_T_async(void **pp, const size_t len, const size_t sizeofT, cudaStream_t stream);
  virtual int deallocate_device_T_async(void **pp, cudaStream_t stream);
  virtual int clear_device_array_T_async(void *data, const size_t ndata, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_HtoD_T(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT);
  virtual int copy_HtoD_T_async(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_DtoH_T(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT);
  virtual int copy_DtoH_T_async(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_DtoD_T(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT);
  virtual int copy_DtoD_T_async(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  /// @}
  /// @brief Functions to get device pointers for atom properties
  /// This functions should be overridden in derived proxy classes that manage actual GPU memory
  /// @{
  virtual float* proxy_atoms_masses_gpu_float() {return nullptr;}
  virtual float* proxy_atoms_charges_gpu_float() {return nullptr;}
  virtual cvm::real* proxy_atoms_masses_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_charges_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_positions_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_total_forces_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_new_colvar_forces_gpu() {return nullptr;}
  /// @}
  /**
   * @brief This function will be called after atom groups are calculated on GPU.
   *
   * This function is useful when additional information is needed to transfer
   * from the proxy. For example, the proxy can copy the lattice vectors in a
   * separate stream, and this function can wait for that stream to complete.
   */
  virtual int wait_for_extra_info_ready();
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)

  /// Pointer to the associated colvarmodule object
  colvarmodule *cvmodule;

  /// Constructor
  colvarproxy();

  /// Destructor
  virtual ~colvarproxy();

  inline std::string const &engine_name() const
  {
    return engine_name_;
  }

  // bool io_available() override;

  /// Request deallocation of the module (currently only implemented by VMD)
  virtual int request_deletion();

  /// Whether deallocation was requested
  inline bool delete_requested() const
  {
    return b_delete_requested;
  }

  /// \brief Reset proxy state, e.g. requested atoms
  virtual int reset();

  /// (Re)initialize the module
  virtual int parse_module_config();

  /// \brief Read a selection of atom IDs from a PDB coordinate file
  /// \param[in] filename name of the file
  /// \param[in,out] atoms array into which atoms will be read from "filename"
  /// \param[in] pdb_field if the file is a PDB and this string is non-empty,
  /// select atoms for which this field is non-zero
  /// \param[in] pdb_field_value if non-zero, select only atoms whose pdb_field equals this
  virtual int load_atoms_pdb(char const *filename, cvm::atom_group &atoms,
                             std::string const &pdb_field, double pdb_field_value);


  /// \brief Load a set of coordinates from a PDB file
  /// \param[in] filename name of the file
  /// \param[in,out] pos array of coordinates to fill; if not empty, the number of its elements must match
  /// the number of entries in "filename"
  /// \param[in] sorted_ids array of sorted internal IDs, used to loop through the file only once
  /// \param[in] pdb_field if non-empty, only atoms for which this field is non-zero will be processed
  /// \param[in] pdb_field_value if non-zero, process only atoms whose pdb_field equals this
  virtual int load_coords_pdb(char const *filename, std::vector<cvm::atom_pos> &pos,
                              std::vector<int> const &sorted_ids, std::string const &pdb_field,
                              double pdb_field_value);

  /// (Re)initialize required member data (called after the module)
  virtual int setup();

  /// Whether the engine allows to fully initialize Colvars immediately
  inline bool engine_ready() const
  {
    return engine_ready_;
  }

  /// Are total forces currently valid? (They would not be right after configuration change)
  inline bool total_forces_valid() const
  {
    return total_forces_valid_;
  }

  /// Mark the total forces as invalid (due to e.g. a configuration change)
  void set_total_forces_invalid();

  /// Mark the total forces as up to date
  void set_total_forces_valid();

  /// Enqueue new configuration text, to be parsed as soon as possible
  void add_config(std::string const &cmd, std::string const &conf);

  /// Update data required by Colvars module (e.g. read atom positions)
  ///
  /// TODO Break up colvarproxy_namd and colvarproxy_lammps function into these
  virtual int update_input();

  /// Update data based on the results of a Colvars call (e.g. send forces)
  virtual int update_output();

  /// Carry out operations needed before next simulation step is run
  int end_of_step();

  /// Print a message to the main log
  virtual void log(std::string const &message);

  /// Print a message to the main log and/or let the host code know about it
  virtual void error(std::string const &message);

  /// Record error message (used by VMD to collect them after a script call)
  void add_error_msg(std::string const &message);

  /// Retrieve accumulated error messages
  std::string const & get_error_msgs();

  /// As the name says
  void clear_error_msgs();

  /// Whether a simulation is running (warn against irrecovarable errors)
  inline bool simulation_running() const
  {
    return b_simulation_running;
  }

  /// Is the current step a repetition of a step just executed?
  /// This is set to true when the step 0 of a new "run" command is being
  /// executed, regardless of whether a state file has been loaded.
  inline bool simulation_continuing() const
  {
    return b_simulation_continuing;
  }

  /// Called at the end of a simulation segment (i.e. "run" command)
  int post_run();

  /// Print a full list of all input atomic arrays for debug purposes
  void print_input_atomic_data();

  /// Print a full list of all applied forces for debug purposes
  void print_output_atomic_data();

  /// Convert a version string "YYYY-MM-DD" into an integer
  int get_version_from_string(char const *version_string);

  /// Get the version number (higher = more recent)
  int version_number() const
  {
    return version_int;
  }

protected:

  /// Next value of lambda to be sent to back-end
  cvm::real cached_alch_lambda;

  /// Whether lambda has been set and needs to be updated in backend
  bool cached_alch_lambda_changed;

  /// Boltzmann constant in internal Colvars units
  cvm::real boltzmann_;

  /// Most up to date target temperature (K units); default to 0.0 if undefined
  cvm::real target_temperature_;

  /// Current integration timestep (engine units); default to 1.0 if undefined
  double timestep_;

  /// Current timestep multiplier, if Colvars is only called once every n MD timesteps
  int time_step_factor_ = 1;

  /// \brief Value of 1 Angstrom in the internal (front-end) Colvars unit for atomic coordinates
  /// * defaults to 0 in the base class; derived proxy classes must set it
  /// * in VMD proxy, can only be changed when no variables are defined
  /// as user-defined values in composite units must be compatible with that system
  cvm::real angstrom_value_;

  /// \brief Value of 1 kcal/mol in the internal Colvars unit for energy
  cvm::real kcal_mol_value_;

  /// Whether the total forces have been requested
  bool total_force_requested;

  /// Use the PBC functions from the Colvars library (as opposed to MD engine)
  bool use_internal_pbc_ = false;

  /// \brief Type of boundary conditions
  ///
  /// Orthogonal and triclinic cells are made available to objects.
  /// For any other conditions (mixed periodicity, triclinic cells in LAMMPS)
  /// minimum-image distances are computed by the host engine regardless.
  enum Boundaries_type {
    boundaries_non_periodic,
    boundaries_pbc_ortho,
    boundaries_pbc_triclinic,
    boundaries_unsupported
  };

  /// Type of boundary conditions
  Boundaries_type boundaries_type;

  /// Bravais lattice vectors
  cvm::rvector unit_cell_x, unit_cell_y, unit_cell_z;

  /// Reciprocal lattice vectors
  cvm::rvector reciprocal_cell_x, reciprocal_cell_y, reciprocal_cell_z;

  /// \brief Array of 0-based integers used to uniquely associate atoms
  /// within the host program
  std::vector<int>          atoms_ids;
  /// \brief Keep track of how many times each atom is used by a separate colvar object
  std::vector<size_t>       atoms_refcount;
  /// \brief Masses of the atoms (allow redefinition during a run, as done e.g. in LAMMPS)
  std::vector<cvm::real, allocator_type<cvm::real>>    atoms_masses;
  /// \brief Charges of the atoms (allow redefinition during a run, as done e.g. in LAMMPS)
  std::vector<cvm::real, allocator_type<cvm::real>>    atoms_charges;
  /// \brief Current three-dimensional positions of the atoms
  std::vector<cvm::rvector, allocator_type<cvm::rvector>> atoms_positions;
  /// \brief Most recent total forces on each atom
  std::vector<cvm::rvector, allocator_type<cvm::rvector>> atoms_total_forces;
  /// \brief Forces applied from colvars, to be communicated to the MD integrator
  std::vector<cvm::rvector, allocator_type<cvm::rvector>> atoms_new_colvar_forces;

  /// Root-mean-square of the applied forces
  cvm::real atoms_rms_applied_force_;

  /// Maximum norm among all applied forces
  cvm::real atoms_max_applied_force_;

  /// ID of the atom with the maximum norm among all applied forces
  int atoms_max_applied_force_id_;

  /// Whether the atom list has been modified internally
  bool modified_atom_list_;

  /// Whether the masses and charges have been updated from the host code
  bool updated_masses_, updated_charges_;

  /// Used by all init_atom() functions: create a slot for an atom not
  /// requested yet; returns the index in the arrays
  int add_atom_slot(int atom_id);

  /// \brief Array of 0-based integers used to uniquely associate atom groups
  /// within the host program
  std::vector<int>          atom_groups_ids;
  /// \brief Keep track of how many times each group is used by a separate cvc
  std::vector<size_t>       atom_groups_refcount;
  /// \brief Total masses of the atom groups
  std::vector<cvm::real>    atom_groups_masses;
  /// \brief Total charges of the atom groups (allow redefinition during a run, as done e.g. in LAMMPS)
  std::vector<cvm::real>    atom_groups_charges;
  /// \brief Current centers of mass of the atom groups
  std::vector<cvm::rvector> atom_groups_coms;
  /// \brief Most recently updated total forces on the com of each group
  std::vector<cvm::rvector> atom_groups_total_forces;
  /// \brief Forces applied from colvars, to be communicated to the MD integrator
  std::vector<cvm::rvector> atom_groups_new_colvar_forces;

  /// Root-mean-square of the applied group forces
  cvm::real atom_groups_rms_applied_force_;

  /// Maximum norm among all applied group forces
  cvm::real atom_groups_max_applied_force_;

  /// Used by all init_atom_group() functions: create a slot for an atom group not requested yet
  int add_atom_group_slot(int atom_group_id);

  /// \brief Array of numeric IDs of volumetric maps
  std::vector<int>          volmaps_ids;

  /// \brief Keep track of how many times each vol map is used by a
  /// separate colvar object
  std::vector<size_t>       volmaps_refcount;

  /// \brief Current values of the vol maps
  std::vector<cvm::real>    volmaps_values;

  /// \brief Forces applied from colvars, to be communicated to the MD
  /// integrator
  std::vector<cvm::real>    volmaps_new_colvar_forces;

  /// Root-mean-square of the the applied forces
  cvm::real volmaps_rms_applied_force_;

  /// Maximum norm among all applied forces
  cvm::real volmaps_max_applied_force_;

  /// Lock state for OpenMP
  omp_lock_t *omp_lock_state;

  /// Whether threaded parallelization should be used (TODO: make this a
  /// cvmodule->deps feature)
  smp_mode_t smp_mode;

  /// MPI communicator containint 1 root proc from each world
  replicas_mpi_comm_t replicas_mpi_comm;

  /// Index (rank) of this replica in the MPI implementation
  int replicas_mpi_rank = 0;

  /// Number of replicas in the MPI implementation
  int replicas_mpi_num = 1;

  /// Pointer to Tcl interpreter object
  Tcl_Interp *tcl_interp_;

  /// Prefix of the input state file to be read next
  std::string input_prefix_str;

  /// Default prefix to be used for all output files (final configuration)
  std::string output_prefix_str;

  /// Prefix of the restart (checkpoint) file to be written next
  std::string restart_output_prefix_str;

  /// How often the simulation engine will write its own restart
  int restart_frequency_engine;

  /// Container of input files/channels indexed by path name
  std::map<std::string, std::istream *> input_streams_;

  /// Object whose reference is returned when read errors occur
  std::istream *input_stream_error_;

  /// Currently open output files/channels
  std::map<std::string, std::ostream *> output_streams_;

  /// Object whose reference is returned when write errors occur
  std::ostream *output_stream_error_;

  /// \brief Whether the proxy supports GPU
  bool support_gpu;

  /// Whether the engine allows to fully initialize Colvars immediately
  bool engine_ready_;

  /// Whether the total forces are currently valid
  bool total_forces_valid_ = false;

  /// Collected error messages
  std::string error_output;

  /// Whether a simulation is running (warn against irrecovarable errors)
  bool b_simulation_running;

  /// Is the current step a repetition of a step just executed?
  /// This is set to true when the step 0 of a new "run" command is being
  /// executed, regardless of whether a state file has been loaded.
  bool b_simulation_continuing;

  /// Whether the entire module should be deallocated by the host engine
  bool b_delete_requested;

  /// Integer representing the version string (allows comparisons)
  int version_int;

  /// Track which features have been acknowledged during the last run
  size_t features_hash;

  /// Name of the simulation engine that the derived proxy object supports
  std::string engine_name_ = "standalone";

  /// Queue of config strings or files to be fed to the module
  void *config_queue_;

};

inline cvm::rvector colvarproxy::position_distance_internal(cvm::atom_pos const &pos1,
                                                            cvm::atom_pos const &pos2) const
{
  if (boundaries_type == boundaries_unsupported) {
    cvm::error_static(cvmodule, "Error: unsupported boundary conditions.\n", COLVARS_INPUT_ERROR);
    return cvm::rvector({0.0, 0.0, 0.0});
  }

  if (boundaries_type == boundaries_non_periodic) {
    return pos2 - pos1;
  }

  // Periodicity flags are hard-coded, because this is the only case supported so far other than
  // the two above
  return position_distance_kernel(pos1, pos2, unit_cell_x, unit_cell_y, unit_cell_z,
                                  reciprocal_cell_x, reciprocal_cell_y, reciprocal_cell_z,
                                  true, true, true);
}


inline cvm::rvector colvarproxy::position_distance_kernel(cvm::atom_pos const &pos1,
                                                                 cvm::atom_pos const &pos2,
                                                                 cvm::rvector const &a,
                                                                 cvm::rvector const &b,
                                                                 cvm::rvector const &c,
                                                                 cvm::rvector const &a_r,
                                                                 cvm::rvector const &b_r,
                                                                 cvm::rvector const &c_r,
                                                                 bool a_p,
                                                                 bool b_p,
                                                                 bool c_p)
{
  cvm::rvector diff = (pos2 - pos1);

  cvm::real const x_shift = std::floor(a_r * diff + 0.5);
  cvm::real const y_shift = std::floor(b_r * diff + 0.5);
  cvm::real const z_shift = std::floor(c_r * diff + 0.5);

  if (a_p) {
    diff.x -= x_shift * a.x + y_shift * b.x + z_shift * c.x;
  }

  if (b_p) {
    diff.y -= x_shift * a.y + y_shift * b.y + z_shift * c.y;
  }

  if (c_p) {
    diff.z -= x_shift * a.z + y_shift * b.z + z_shift * c.z;
  }

  return diff;
}

#endif
