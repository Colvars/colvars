// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_NAMD_H
#define COLVARPROXY_NAMD_H

#include <memory>

#include "colvarproxy_namd_version.h"

// For NAMD_UNIFIED_REDUCTION and AtomIDList
#include "NamdTypes.h"
// For CMK_SMP && USE_CKLOOP
#include "Node.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"


class Controller;
class GlobalMasterColvars;
class GridforceFullMainGrid;
class Random;
class SimParameters;
class Molecule;
class ScriptTcl;
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
class CudaGlobalMasterColvars;
#endif


/// Communication between colvars and NAMD (implementation of \link colvarproxy \endlink)
class colvarproxy_namd : public colvarproxy {
public:

  colvarproxy_namd();
  ~colvarproxy_namd();

  /// Tell the proxy that it will be using the GlobalMaster interface
  void set_gm_object(GlobalMasterColvars *gm);

  /// Initialize proxy data based on the chosen interface
  int init(); // no override

  /// Initialize the module
  int init_module(); // no override

  void init_tcl_pointers() override;
  int setup() override;
  int reset() override;

  /// Get the target temperature from the NAMD thermostats supported so far
  int update_target_temperature();

  /// Create map from NAMD atom indices to colvarproxy array indices (used by GlobalMaster)
  void init_gm_atoms_map();

  /// Request the given atom through the GlobalMaster object
  /// \param aid Atom ID to request
  /// \param index Index in the colvarproxy array
  void request_gm_atom_by_id(int aid, int index);

  /// Update map to reflect other requested atoms, including other GlobalMaster objects
  int update_gm_atoms_map(AtomIDList::const_iterator begin, AtomIDList::const_iterator end);

  /// Create and zero out data buffers for atoms requested through GlobalMaster
  int setup_gm_atom_buffers();

  /// Create and zero out data buffers for atom groups requested through GlobalMaster
  int setup_gm_atom_group_buffers();

  /// Create and zero out data buffers for grid objects requested through GlobalMaster
  int setup_gm_volmap_buffers();

  /// Read the current coordinates and total forces
  void read_gm_atom_buffers();

  /// Send Colvars forces to GlobalMaster
  void send_gm_atom_forces();

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
  /// Initialize on the CUDA worker, using engine objects supplied by the client.
  int initialize_from_cudagm(CudaGlobalMasterColvars *client,
                             std::vector<std::string> const &arguments,
                             int deviceID, cudaStream_t stream,
                             SimParameters const *parameters,
                             Molecule const *molecule, ScriptTcl *script);
  bool atomsChanged() const { return modified_atom_list(); }
  void reallocate();
  void onBuffersUpdated();

  cudaStream_t getStream() const { return mStream; }
  double getEnergy() const { return mBiasEnergy; }
  double *getPositions() const { return d_mPositions; }
  double *getAppliedForces() const { return d_mAppliedForces; }
  double *getTotalForces() const { return d_mTotalForces; }
  float *getMasses() const { return d_mMass; }
  float *getCharges() const { return d_mCharges; }
  double *getLattice() const { return d_mLattice; }

  cudaStream_t get_default_stream() override { return mStream; }
  float *proxy_atoms_masses_gpu_float() override { return d_mMass; }
  float *proxy_atoms_charges_gpu_float() override { return d_mCharges; }
  cvm::real *proxy_atoms_positions_gpu() override { return d_mPositions; }
  cvm::real *proxy_atoms_total_forces_gpu() override { return d_mTotalForces; }
  cvm::real *proxy_atoms_new_colvar_forces_gpu() override { return d_mAppliedForces; }
  cvm::system_boundary_conditions get_system_boundaries() override;

  friend class CudaGlobalMasterColvars;
#endif

protected:

  bool has_cudagm_client() const
  {
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
    return mClient != nullptr;
#else
    return false;
#endif
  }
  /// CUDA workers must not resolve engine data through Node::Object().
  Molecule *get_molecule() const;

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
  int allocateDeviceArrays();
  int deallocateDeviceArrays();
  int allocateDeviceTransposeArrays();
  int deallocateDeviceTransposeArrays();
  void read_cudagm_buffers();
  void send_cudagm_forces();
  void set_lattice();

  CudaGlobalMasterColvars *mClient = nullptr;
  Molecule *mMolecule = nullptr;
  ScriptTcl *mScriptTcl = nullptr;
  int m_device_id = -1;
  cudaStream_t mStream = nullptr;
  double mBiasEnergy = 0.0;
  double *d_mPositions = nullptr;
  double *d_mAppliedForces = nullptr;
  double *d_mTotalForces = nullptr;
  float *d_mMass = nullptr;
  float *d_mCharges = nullptr;
  double *d_mLattice = nullptr;
  double *h_mLattice = nullptr;
  cvm::rvector *d_trans_mPositions = nullptr;
  cvm::rvector *d_trans_mAppliedForces = nullptr;
  cvm::rvector *d_trans_mTotalForces = nullptr;
  cvm::real *d_trans_mMass = nullptr;
  cvm::real *d_trans_mCharges = nullptr;
  size_t allocated_atoms = 0;
  bool lattice_copy_pending = false;
#endif

  /// Pointer to the parent GlobalMaster object
  GlobalMasterColvars *globalmaster = nullptr;

  /// Map from NAMD atom indices to colvarproxy array indices (used by GlobalMaster)
  std::vector<int> gm_atoms_map;

  /// Pointer to the NAMD simulation input object
  SimParameters *simparams = nullptr;

  /// Pointer to Controller object
  Controller const *controller = nullptr;

  /// NAMD-style PRNG object
  std::unique_ptr<Random> random;

  /// Use to distinguish between "run 0" and actual runs
  bool first_timestep = true;

  /// Current NAMD simulation step (promoted from int)
  cvm::step_number NAMD_step = 0L;

  /// Previous NAMD simulation step; used to test if the simulation is advancing
  cvm::step_number previous_NAMD_step = 0L;

  /// Used to submit restraint energy as MISC
#if !defined (NAMD_UNIFIED_REDUCTION)
  SubmitReduction *reduction = nullptr;
#endif
#if defined(NODEGROUP_FORCE_REGISTER) && !defined(NAMD_UNIFIED_REDUCTION)
  NodeReduction *nodeReduction = nullptr;
#endif

  /// Accelerated MD reweighting factor
  bool accelMDOn = false;
  cvm::real amd_weight_factor = 1.0;
  void update_accelMD_info();

public:

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

  smp_mode_t get_preferred_smp_mode() const override;
  std::vector<smp_mode_t> get_available_smp_modes() const override;
  int set_smp_mode(smp_mode_t mode) override;

#if CMK_SMP && USE_CKLOOP
  int smp_loop(int n_items, std::function<int (int)> const &worker) override;

  int smp_biases_loop() override;

  int smp_biases_script_loop() override;

  friend void calc_colvars_items_smp(int first, int last, void *result, int paramNum, void *param);
  friend void calc_cv_biases_smp(int first, int last, void *result, int paramNum, void *param);
  friend void calc_cv_scripted_forces(int paramNum, void *param);

  int smp_thread_id()
  {
    return has_cudagm_client() ? 0 : CkMyRank();
  }

  int smp_num_threads()
  {
    return has_cudagm_client() ? 1 : CkMyNodeSize();
  }

protected:

  CmiNodeLock charm_lock_state;
  bool charm_lock_initialized = false;

public:

  int smp_lock()
  {
    if (charm_lock_initialized) CmiLock(charm_lock_state);
    return COLVARS_OK;
  }

  int smp_trylock()
  {
    const int ret = charm_lock_initialized ? CmiTryLock(charm_lock_state) : 0;
    if (ret == 0) return COLVARS_OK;
    else return COLVARS_ERROR;
  }

  int smp_unlock()
  {
    if (charm_lock_initialized) CmiUnlock(charm_lock_state);
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

  enum class e_pdb_field {
    none,
    occ,
    beta,
    x,
    y,
    z,
    ntot
  };

  e_pdb_field pdb_field_str2enum(std::string const &pdb_field_str);

  int load_atoms_pdb(char const *filename,
                     cvm::atom_group &atoms,
                     std::string const &pdb_field,
                     double pdb_field_value) override;

  int load_coords_pdb(char const *filename,
                      std::vector<cvm::atom_pos> &pos,
                      const std::vector<int> &indices,
                      std::string const &pdb_field,
                      double const pdb_field_value) override;


  int check_scalable_group_coms() override;

  int init_atom_group(std::vector<int> const &atoms_ids) override;
  void clear_atom_group(int index) override;

  int update_group_properties(int index);

  int check_volmaps_available() override;

  int check_engine_volmaps_available() override;


  /// Select a MGridForces map for computation by NAMD
  int request_engine_volmap_by_id(int volmap_id) override;

  /// Select a MGridForces map for computation by NAMD
  int request_engine_volmap_by_name(std::string const &volmap_name) override;

  /// Add map to GlobalMaster client (if not already in them)
  void request_globalmaster_volmap(int volmap_id);

  /// Select a MGridForces map for internal computation (frontend)
  int init_internal_volmap_by_id(int volmap_id) override;

  /// Select a MGridForces map for internal computation (frontend)
  int init_internal_volmap_by_name(std::string const &volmap_name) override;

  /// Load a map internally independent from MGridForces
  int load_internal_volmap_from_file(std::string const &filename) override;

  int clear_volmap(int index) override;

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
