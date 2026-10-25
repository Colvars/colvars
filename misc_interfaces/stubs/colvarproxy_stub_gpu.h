// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_STUBS_GPU_H
#define COLVARPROXY_STUBS_GPU_H


#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
class colvarproxy_stub_gpu : public colvarproxy {
public:
  colvarproxy_stub_gpu();
  ~colvarproxy_stub_gpu() override;
  int setup() override;
  void request_total_force(bool yesno) override;
  bool total_forces_enabled() const override;
  bool total_forces_same_step() const override;
  void log(std::string const &message) override;
  void error(std::string const &message) override;
  int set_unit_system(std::string const &units_in, bool check_only) override;
  int init_atom(int atom_number) override;
  int check_atom_id(int atom_number) override;
  void clear_atom(int index) override;
  int read_frame_xyz(const char *filename, const bool write_force_file = false);
  void reallocate() {
    deallocateDeviceArrays();
    allocateDeviceArrays();
  }
  // float* proxy_atoms_masses_gpu_float() override {return d_mMass;}
  // float* proxy_atoms_charges_gpu_float() override {return d_mCharges;}
  cvm::real* proxy_atoms_positions_gpu() override {return d_mPositions;}
  cvm::real* proxy_atoms_total_forces_gpu() override {return d_mTotalForces;}
  cvm::real* proxy_atoms_new_colvar_forces_gpu() override {return d_mAppliedForces;}
  cudaStream_t get_default_stream() override {return stream;}
  smp_mode_t get_preferred_smp_mode() const override {
    return smp_mode_t::gpu;
  }
  std::vector<smp_mode_t> get_available_smp_modes() const override {
    std::vector<colvarproxy_smp::smp_mode_t> available_modes{
      smp_mode_t::gpu
    };
    return available_modes;
  }
  int set_smp_mode(smp_mode_t mode) override {
    if (mode == smp_mode_t::gpu) {
      smp_mode = mode;
      support_gpu = true;
    }
    return COLVARS_NOT_IMPLEMENTED;
  }
  void init_cvm() {
    cvmodule = new colvarmodule(this);
    cvmodule->log("Using minimal CUDA testing interface.\n");

    cvmodule->cv_traj_freq = 0; // I/O will be handled explicitly
    cvmodule->restart_out_freq = 0;
    cvm::rotation::monitor_crossings = false; // Avoid unnecessary error messages

    cvmodule->setup_input();
    cvmodule->setup_output();

    colvarproxy_stub_gpu::setup();
  }
private:
  void allocateDeviceArrays() {
    const int numAtoms = atoms_ids.size();
    allocate_device(&d_mPositions, 3*numAtoms);
    allocate_device(&d_mAppliedForces, 3*numAtoms);
    allocate_device(&d_mTotalForces, 3*numAtoms);
  }
  void deallocateDeviceArrays() {
    deallocate_device(&d_mPositions);
    deallocate_device(&d_mAppliedForces);
    deallocate_device(&d_mTotalForces);
  }
  cudaStream_t stream;
  double* d_mPositions;
  double* d_mAppliedForces;
  double* d_mTotalForces;
  bool mAtomsChanged;
  // These buffers are for read_frame_xyz
  std::vector<cvm::rvector> m_positions;
  colvarproxy_atoms::atom_buffer_real_t m_positions_soa;
  colvarproxy_atoms::atom_buffer_real_t m_applied_forces_soa;
  std::vector<cvm::rvector> m_forces_by_id;
};

#endif
#endif
