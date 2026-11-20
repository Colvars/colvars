#include "CudaGlobalMasterClient.h"
#include "colvar_gpu_support.h"
#include "colvarproxy_cudaglobalmaster.h"
#include "colvarproxy_cudaglobalmaster_kernel.h"
#include "colvarproxy_namd.h"
#include "Molecule.h"
#include "CudaUtils.h"
#include "Random.h"
#include "colvarparse.h"
#include "colvaratoms.h"
#include "ScriptTcl.h"
#include "colvarscript.h"

#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
#include <nvtx3/nvToolsExt.h>
#include <nvtx3/nvToolsExtCuda.h>
#endif // CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING

#if defined (NAMD_CUDA) || defined (NAMD_HIP)

#if defined (__linux__) || defined (__APPLE__)
extern "C" {
  CudaGlobalMasterColvars* allocator() {
    return new CudaGlobalMasterColvars();
  }
  void deleter(CudaGlobalMasterColvars* ptr) {
    if (ptr != nullptr) {
      delete ptr;
    }
  }
}
#endif

#ifdef WIN32
extern "C" {
  __declspec (dllexport) CudaGlobalMasterColvars* allocator() {
    return new CudaGlobalMasterColvars();
  }
  __declspec (dllexport) void deleter(CudaGlobalMasterColvars* ptr) {
    if (ptr != nullptr) {
      delete ptr;
    }
  }
}
#endif

class colvarproxy_impl: public colvarproxy_namd {
public:
  colvarproxy_impl();
  virtual ~colvarproxy_impl();
  void add_energy(cvm::real energy) override { mBiasEnergy += energy; }
  void clear_atom(int index) override;
  int setup() override;
  /// Volumetric maps are not yet supported by this interface
  int check_volmaps_available() override { return COLVARS_NOT_IMPLEMENTED; }
  int check_replicas_enabled() override;
  int replica_index() override;
  int num_replicas() override;
  void replica_comm_barrier() override;
  int replica_comm_recv(char* msg_data, int buf_len, int src_rep) override;
  int replica_comm_send(char* msg_data, int msg_len, int dest_rep) override;
  void initialize_from_cudagm(
    CudaGlobalMasterColvars *client,
    SimParameters *s,
    Molecule *m,
    const std::vector<std::string> &arguments,
    const int deviceID, cudaStream_t stream);
  const bool atomsChanged() const { return modified_atom_list(); }
  void calculate() override;
  void copy_cudagm_data_from_device();
  void copy_cudagm_data_to_device();
  void onBuffersUpdated();
  // FIXME: This interface does not run Colvars on the same thread as ScriptTcl,
  //        so it is buggy to use the same interpreter as NAMD.
  void init_tcl_pointers() override {
    colvarproxy::init_tcl_pointers(); // Create dedicated interpreter
  }
  int reset() override;
  void reallocate() {
    deallocateDeviceArrays();
    allocateDeviceArrays();
    deallocateDeviceTransposeArrays();
    allocateDeviceTransposeArrays();
    cvmodule->proxy_gpu_buffers_reallocated_done();
  }
  smp_mode_t get_preferred_smp_mode() const override {
    return smp_mode_t::none;
  }
  std::vector<smp_mode_t> get_available_smp_modes() const override {
    std::vector<colvarproxy_smp::smp_mode_t> available_modes{
      smp_mode_t::none,
      smp_mode_t::gpu
    };
    return available_modes;
  }
  int set_smp_mode(smp_mode_t mode) override {
    if ((mode == smp_mode_t::none) || (mode == smp_mode_t::gpu)) {
      smp_mode = mode;
      if (mode == smp_mode_t::gpu) {
        support_gpu = true;
#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
        mEventAttrib.message.ascii = "Colvars GPU";
#endif
      } else {
#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
        mEventAttrib.message.ascii = "Colvars CPU";
#endif
        support_gpu = false;
      }
      return COLVARS_OK;
    }
    return COLVARS_NOT_IMPLEMENTED;
  }
  float* proxy_atoms_masses_gpu_float() override {return d_mMass;}
  float* proxy_atoms_charges_gpu_float() override {return d_mCharges;}
  cvm::real* proxy_atoms_positions_gpu() override {return d_mPositions;}
  cvm::real* proxy_atoms_total_forces_gpu() override {return d_mTotalForces;}
  cvm::real* proxy_atoms_new_colvar_forces_gpu() override {return d_mAppliedForces;}
  cudaStream_t get_default_stream() override {return mStream;}
  void set_lattice();
  cvm::system_boundary_conditions get_system_boundaries()  override {
    if (mClient->requestUpdateLattice()) {
      if (has_gpu_support()) {
        cudaCheck(cudaEventSynchronize(get_event(colvarproxy_gpu::event_type::update_lattice)));
      }
      set_lattice();
    }
    return boundaries_;
  };
  friend class CudaGlobalMasterColvars;
private:
  void allocateDeviceArrays();
  void deallocateDeviceArrays();
  void allocateDeviceTransposeArrays();
  void deallocateDeviceTransposeArrays();
  double* d_mPositions;
  double* d_mAppliedForces;
  double* d_mTotalForces;
  double* d_mLattice;
  float*  d_mMass;
  float*  d_mCharges;
  // For transpose
  cvm::rvector* d_trans_mPositions;
  cvm::rvector* d_trans_mAppliedForces;
  cvm::rvector* d_trans_mTotalForces;
  cvm::real*    d_trans_mMass;
  cvm::real*    d_trans_mCharges;
  double* h_mLattice;
  double mBiasEnergy;
  cudaStream_t mStream;
  std::vector<std::string> mConfigFiles;
  CudaGlobalMasterColvars* mClient;
  int m_device_id;
#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
  nvtxEventAttributes_t mEventAttrib;
#endif // CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
};


colvarproxy_impl::colvarproxy_impl(): colvarproxy_namd(),
  d_mPositions(nullptr), d_mAppliedForces(nullptr),
  d_mTotalForces(nullptr), d_mLattice(nullptr),
  d_mMass(nullptr), d_mCharges(nullptr),
  d_trans_mPositions(nullptr),
  d_trans_mAppliedForces(nullptr),
  d_trans_mTotalForces(nullptr),
  d_trans_mMass(nullptr),
  d_trans_mCharges(nullptr),
  h_mLattice(nullptr),
  mBiasEnergy(0), mStream(0) {
  // colvars = nullptr;
  cvmodule = nullptr;
#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
  mEventAttrib.version = NVTX_VERSION;
  mEventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  mEventAttrib.colorType = NVTX_COLOR_ARGB;
  mEventAttrib.color = 0xFF880000;
  mEventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
  mEventAttrib.message.ascii = "Colvars CPU";
#endif // CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
}

colvarproxy_impl::~colvarproxy_impl() {
  deallocateDeviceArrays();
  deallocateDeviceTransposeArrays();
}

int colvarproxy_impl::reset() {
  int error_code = COLVARS_OK;
  deallocateDeviceArrays();
  deallocateDeviceTransposeArrays();
  modified_atom_list_ = true;
  error_code |= colvarproxy_namd::reset();
  return error_code;
}

int colvarproxy_impl::setup() {
  colvarproxy::parse_module_config();
  int error_code = colvarproxy::setup();
  if (cvmodule->size() == 0) {
    // Module is empty, nothing to do
    return COLVARS_OK;
  }
  if (simparams->wrapAll) {
    log("Warning: enabling wrapAll can lead to inconsistent results "
        "for Colvars calculations: please disable wrapAll, "
        "as is the default option in NAMD.\n");
  }
  set_integration_timestep(simparams->dt);
  log("updating target temperature (T = "+
      cvm::to_str(target_temperature())+" K).\n");
  error_code |= cvmodule->update_engine_parameters();
  error_code |= cvmodule->setup_input();
  error_code |= cvmodule->setup_output();
  if (has_gpu_support()) {
    cvmodule->log("Warning: the GPU support in Colvars is still experimental!\n");
  }
  return error_code;
}

void colvarproxy_impl::initialize_from_cudagm(
  CudaGlobalMasterColvars* client,
  SimParameters* s,
  Molecule* m,
  const std::vector<std::string>& arguments,
  const int deviceID, cudaStream_t stream) {
  mClient = client;
  simparams = s;
  molecule = m;
  random.reset(new Random(simparams->randomSeed));
  mStream = stream;
  m_device_id = deviceID;
  int savedDevice;
  cudaCheck(cudaGetDevice(&savedDevice));
  cudaCheck(cudaSetDevice(m_device_id));
  allocate_device<double>(&d_mLattice, sizeof(double)*4*3);
  allocate_host<double>(&h_mLattice, sizeof(double)*4*3);
  cudaCheck(cudaSetDevice(savedDevice));
  if (arguments.size() < 3) {
    const std::string error = "Wrong number of arguments of CudaGlobalMasterColvars.";
    NAMD_die(error.c_str());
  }
  mConfigFiles.clear();
  mConfigFiles.insert(mConfigFiles.end(), arguments.begin()+2, arguments.end());
  // const int64_t step = mClient->getStep();
  engine_name_ = "NAMD_CUDAGLOBALMASTER";
  // if (colvars != nullptr) delete colvars;
  if (cvmodule == nullptr) {
    // std::cerr << "Creating cvmodule" << std::endl;
    // printf("Creating cvmodule\n");
    cvmodule = new colvarmodule(this);
  }
  cvmodule->log("Using " + engine_name_ + " interface, version "+
           cvm::to_str(0)+".\n");
  cvmodule->cite_feature("NAMD engine");
  cvmodule->cite_feature("Colvars-NAMD interface");
  for (auto it = mConfigFiles.begin(); it != mConfigFiles.end(); ++it) {
    add_config("configfile", *it);
  }
  update_target_temperature();
  setup();
  if (simparams->firstTimestep != 0) {
    cvmodule->set_initial_step(static_cast<cvm::step_number>(simparams->firstTimestep));
  }
  colvarproxy_io::set_output_prefix(std::string(simparams->outputFilename));
  colvarproxy_io::set_restart_output_prefix(std::string(simparams->restartFilename));
  colvarproxy_io::set_default_restart_frequency(simparams->restartFrequency);
}

void colvarproxy_impl::clear_atom(int index) {
  colvarproxy_namd::clear_atom(index);
  modified_atom_list_ = true;
}

void colvarproxy_impl::allocateDeviceArrays() {
  const int numAtoms = atoms_ids.size();
  int savedDevice;
  cudaCheck(cudaGetDevice(&savedDevice));
  cudaCheck(cudaSetDevice(m_device_id));
  allocate_device<double>(&d_mPositions, 3*numAtoms);
  allocate_device<double>(&d_mAppliedForces, 3*numAtoms);
  if (mClient->requestedTotalForcesAtomsChanged()) {
    allocate_device<double>(&d_mTotalForces, 3*numAtoms);
  }
  allocate_device<float>(&d_mMass, numAtoms);
  allocate_device<float>(&d_mCharges, numAtoms);
  cudaCheck(cudaSetDevice(savedDevice));
}

void colvarproxy_impl::deallocateDeviceArrays() {
  int savedDevice;
  cudaCheck(cudaGetDevice(&savedDevice));
  cudaCheck(cudaSetDevice(m_device_id));
  deallocate_device<double>(&d_mPositions);
  deallocate_device<double>(&d_mAppliedForces);
  if (mClient->requestedTotalForcesAtomsChanged()) {
    deallocate_device<double>(&d_mTotalForces);
  }
  deallocate_device<float>(&d_mMass);
  deallocate_device<float>(&d_mCharges);
  cudaCheck(cudaSetDevice(savedDevice));
}

void colvarproxy_impl::allocateDeviceTransposeArrays() {
  const int numAtoms = atoms_ids.size();
  allocate_device<cvm::rvector>(&d_trans_mPositions, numAtoms);
  allocate_device<cvm::rvector>(&d_trans_mAppliedForces, numAtoms);
  if (mClient->requestedTotalForcesAtomsChanged()) {
    allocate_device<cvm::rvector>(&d_trans_mTotalForces, numAtoms);
  }
  allocate_device<cvm::real>(&d_trans_mMass, numAtoms);
  allocate_device<cvm::real>(&d_trans_mCharges, numAtoms);
}

void colvarproxy_impl::deallocateDeviceTransposeArrays() {
  deallocate_device<cvm::rvector>(&d_trans_mPositions);
  deallocate_device<cvm::rvector>(&d_trans_mAppliedForces);
  if (mClient->requestedTotalForcesAtomsChanged()) {
    deallocate_device<cvm::rvector>(&d_trans_mTotalForces);
  }
  deallocate_device<cvm::real>(&d_trans_mMass);
  deallocate_device<cvm::real>(&d_trans_mCharges);
}


int colvarproxy_impl::check_replicas_enabled() {
  return mClient->replica_enabled();
}

int colvarproxy_impl::replica_index() {
  return mClient->replica_index();
}

int colvarproxy_impl::num_replicas() {
  return mClient->num_replicas();
}

void colvarproxy_impl::replica_comm_barrier() {
  return mClient->replica_comm_barrier();
}

int colvarproxy_impl::replica_comm_recv(char* msg_data, int buf_len, int src_rep) {
  return mClient->replica_comm_recv(msg_data, buf_len, src_rep);
}

int colvarproxy_impl::replica_comm_send(char* msg_data, int msg_len, int dest_rep) {
  return mClient->replica_comm_send(msg_data, msg_len, dest_rep);
}

void colvarproxy_impl::onBuffersUpdated() {
  const size_t numAtoms = atoms_ids.size();
  int savedDevice;
  cudaCheck(cudaGetDevice(&savedDevice));
  cudaCheck(cudaSetDevice(m_device_id));
  if (!has_gpu_support()) {
    // Clear the previous applied forces
    auto &colvars_applied_force = *(modify_atom_applied_forces());
    // TODO: Why do I need to clean the applied forces manually?
    std::fill(colvars_applied_force.begin(),
              colvars_applied_force.end(),
              cvm::rvector(0, 0, 0));
  } else {
    cudaCheck(cudaMemsetAsync(proxy_atoms_new_colvar_forces_gpu(), 0, 3 * numAtoms * sizeof(double), mStream));
  }
  // Clear the previous bias energy
  mBiasEnergy = 0;
  bool sync = false;
  if (numAtoms > 0) {
    if (!has_gpu_support()) {
      transpose_to_host_rvector(d_mPositions, d_trans_mPositions, numAtoms, mStream);
      if (mClient->requestUpdateAtomTotalForces()) {
        transpose_to_host_rvector(d_mTotalForces, d_trans_mTotalForces, numAtoms, mStream);
      }
      sync = true;
    }
    if (mClient->requestUpdateMasses()) {
      copy_float_to_host_double(d_mMass, d_trans_mMass, numAtoms, mStream);
      sync = true;
    }
    if (mClient->requestUpdateCharges()) {
      copy_float_to_host_double(d_mCharges, d_trans_mCharges, numAtoms, mStream);
      sync = true;
    }
  }
  if (sync) {
    // Check if CUDA kernels are successfully executed
    cudaCheck(cudaPeekAtLastError());
    // Ensure the above kernels are done before the NB forces
    cudaCheck(cudaStreamSynchronize(mStream));
  }
  // Restore the GPU device
  cudaCheck(cudaSetDevice(savedDevice));
}


void colvarproxy_impl::copy_cudagm_data_from_device() {
  const size_t & numAtoms = atoms_ids.size();
  cudaCheck(cudaSetDevice(m_device_id));
  if (mClient->requestUpdateLattice()) {
    ::copy_DtoH(d_mLattice, h_mLattice, 3*4, mStream);
    if (has_gpu_support()) {
      cudaCheck(cudaEventRecord(get_event(colvarproxy_gpu::event_type::update_lattice), mStream));
    }
  }
  // The following memcpy operations are supposed to be overlapped with the NB kernel
  if (numAtoms > 0) {
#if CUDAGM_VERSION >= 3
    if (mClient->requestUpdateAtomTotalForces()) {
      if (mClient->isStartupStep()) {
        set_total_forces_invalid();
      } else {
        set_total_forces_valid();
      }
    }
#endif
    if (!has_gpu_support()) {
      // Transform the arrays for Colvars
      auto &colvars_pos = *(modify_atom_positions());
      ::copy_DtoH(d_trans_mPositions, colvars_pos.data(), numAtoms, mStream);
      if (mClient->requestUpdateAtomTotalForces()) {
        auto &colvars_total_force = *(modify_atom_total_forces());
        ::copy_DtoH(d_trans_mTotalForces, colvars_total_force.data(), numAtoms, mStream);
      }
    }
    if (mClient->requestUpdateMasses()) {
      auto &colvars_mass = *(modify_atom_masses());
      ::copy_DtoH(d_trans_mMass, colvars_mass.data(), numAtoms, mStream);
    }
    if (mClient->requestUpdateCharges()) {
      auto &colvars_charge  = *(modify_atom_charges());
      ::copy_DtoH(d_trans_mCharges, colvars_charge.data(), numAtoms, mStream);
    }
    if (has_gpu_support()) {
      cudaCheck(cudaEventRecord(get_event(colvarproxy_gpu::event_type::copy_atoms), mStream));
    }
  }
  // NOTE: I think the implementation in Colvars will syncrhonize the stream before
  // calculating CVCs anyway, so I can skip it here.
  if (!has_gpu_support()) {
    // Synchronize the stream to make sure the host buffers are ready
    cudaCheck(cudaStreamSynchronize(mStream));
  }
}

void colvarproxy_impl::calculate() {
  NAMD_step = mClient->getStep();
  int savedDevice;
  cudaCheck(cudaGetDevice(&savedDevice));
  if (first_timestep) {
    // setup();
    update_target_temperature();
    cvmodule->update_engine_parameters();
    cvmodule->setup_input();
    cvmodule->setup_output();
    first_timestep = false;
  } else {
    if ( NAMD_step - previous_NAMD_step == 1 ) {
      cvmodule->it++;
      b_simulation_continuing = false;
    } else {
      // Cases covered by this condition:
      // - run 0
      // - beginning of a new run statement
      // The internal counter is not incremented, and the objects are made
      // aware of this via the following flag
      b_simulation_continuing = true;
      // Update NAMD output and restart prefixes
      colvarproxy_io::set_output_prefix(std::string(simparams->outputFilename));
      colvarproxy_io::set_restart_output_prefix(std::string(simparams->restartFilename));
      colvarproxy_io::set_default_restart_frequency(simparams->restartFrequency);
      cvmodule->setup_output();
    }
  }
  previous_NAMD_step = NAMD_step;
  copy_cudagm_data_from_device();

  // Run Colvars
#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
  nvtxRangePushEx(&mEventAttrib);
#endif // CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
  if (cvm::debug()) {
    print_input_atomic_data();
  }
  if (cvmodule->calc() != COLVARS_OK) {
    cvmodule->error("Error in the collective variables module.\n", COLVARS_ERROR);
  }
  if (cvm::debug()) {
    print_output_atomic_data();
  }
#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
  nvtxRangePop();
#endif // CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING

  copy_cudagm_data_to_device();

  // NAMD does not destruct GlobalMaster objects, so we must remember
  // to write all output files at the end of a run
  if (NAMD_step == simparams->N) {
#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
  nvtxRangePushEx(&mEventAttrib);
#endif // CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
    post_run();
#ifdef CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
  nvtxRangePop();
#endif // CUDAGLOBALMASTERCOLVARS_CUDA_PROFILING
  }
  // Restore the GPU device
  cudaCheck(cudaSetDevice(savedDevice));
  // Reset the modified atom list flag
  reset_modified_atom_list();
}

void colvarproxy_impl::copy_cudagm_data_to_device()
{
  const size_t & numAtoms = atoms_ids.size();
  // Update applied forces
  if (!has_gpu_support()) {
    auto &colvars_applied_force = *(modify_atom_applied_forces());
    if (numAtoms > 0) {
      ::copy_HtoD(colvars_applied_force.data(), d_trans_mAppliedForces, numAtoms, mStream);
      transpose_from_host_rvector(
        d_mAppliedForces, d_trans_mAppliedForces,
        numAtoms, mStream);
    }
  }
  // NOTE: I think I can skip the syncrhonization here because this client
  //       share the same stream as the CudaGlobalMasterServer object
  // cudaCheck(cudaStreamSynchronize(mStream));
}

void colvarproxy_impl::set_lattice() {
  const cvm::rvector a1(h_mLattice[0], h_mLattice[1], h_mLattice[2]);
  const cvm::rvector a2(h_mLattice[3], h_mLattice[4], h_mLattice[5]);
  const cvm::rvector a3(h_mLattice[6], h_mLattice[7], h_mLattice[8]);
  const int p1 = ( a1.norm2() ? 1 : 0 );
  const int p2 = ( a2.norm2() ? 1 : 0 );
  const int p3 = ( a3.norm2() ? 1 : 0 );
  boundaries_.set_boundaries(p1, p2, p3, a1, a2, a3);
}

CudaGlobalMasterColvars::CudaGlobalMasterColvars():
  CudaGlobalMasterClient()
{
  if (CudaGlobalMasterClient::getSimParameters()->colvarsOn) {
    NAMD_die("This plugin is incompatible with the Colvars bundled with NAMD.");
  }
  mImpl = std::make_unique<colvarproxy_impl>();
}

CudaGlobalMasterColvars::~CudaGlobalMasterColvars() {}

void CudaGlobalMasterColvars::initialize(
  const std::vector<std::string>& arguments,
  int deviceID, cudaStream_t stream) {
  CudaGlobalMasterClient::initialize(arguments, deviceID, stream);
  mImpl->initialize_from_cudagm(this,
                                CudaGlobalMasterClient::getSimParameters(),
                                CudaGlobalMasterClient::getMolecule(),
                                arguments, deviceID, stream);
}

bool CudaGlobalMasterColvars::requestedAtomsChanged()  {
  return mImpl->atomsChanged();
}

const std::vector<AtomID>& CudaGlobalMasterColvars::getRequestedTotalForcesAtoms() const {
  if (mImpl->total_forces_enabled()) {
    return *(mImpl->get_atom_ids());
  } else {
    return mEmpty;
  }
}

bool CudaGlobalMasterColvars::requestedTotalForcesAtomsChanged() {
  if (mImpl->total_forces_enabled()) {
    return mImpl->atomsChanged();
  } else {
    return false;
  }
}

bool CudaGlobalMasterColvars::requestUpdateMasses() {
  return mImpl->atomsChanged();
}

bool CudaGlobalMasterColvars::requestUpdateCharges() {
  return mImpl->atomsChanged();
}

#if CUDAGM_VERSION >= 3
void CudaGlobalMasterColvars::setStep(int64_t step, int startup, int doMigration) {
  CudaGlobalMasterClient::setStep(step, startup, doMigration);
#else
void CudaGlobalMasterColvars::setStep(int64_t step) {
  CudaGlobalMasterClient::setStep(step);
#endif
  if (mImpl->atomsChanged()) {
    mImpl->reallocate();
  }
}

void CudaGlobalMasterColvars::calculate() {
  mImpl->calculate();
}

cudaStream_t CudaGlobalMasterColvars::getStream() {
  return mImpl->mStream;
}

bool CudaGlobalMasterColvars::requestUpdateAtomTotalForces() {
  return mImpl->total_forces_enabled();
}

double CudaGlobalMasterColvars::getEnergy() const {
  return mImpl->mBiasEnergy;
}

double* CudaGlobalMasterColvars::getAppliedForces() const {
  return mImpl->d_mAppliedForces;
}

double* CudaGlobalMasterColvars::getPositions() {
  return mImpl->d_mPositions;
}

float* CudaGlobalMasterColvars::getMasses() {
  return mImpl->d_mMass;
}

float* CudaGlobalMasterColvars::getCharges() {
  return mImpl->d_mCharges;
}

double* CudaGlobalMasterColvars::getTotalForces() {
  return mImpl->d_mTotalForces;
}

double* CudaGlobalMasterColvars::getLattice() {
  return mImpl->d_mLattice;
}

const std::vector<AtomID>& CudaGlobalMasterColvars::getRequestedAtoms() const {
  return *(mImpl->get_atom_ids());
}

void CudaGlobalMasterColvars::onBuffersUpdated() {
  mImpl->onBuffersUpdated();
}

int CudaGlobalMasterColvars::updateFromTCLCommand(const std::vector<std::string>& arguments) {
  // Prepare the arguments
  const int objc = arguments.size();
  unsigned char** objv = new unsigned char*[objc];
  for (int i = 0; i < objc; ++i) {
    const int len = std::strlen(arguments[i].c_str());
    objv[i] = new unsigned char[len + 1];
    std::strncpy(reinterpret_cast<char*>(objv[i]),
                 arguments[i].c_str(), len + 1);
    objv[i][len] = '\0';
  }
  // Call Colvars scripting interface
  int error_code = TCL_ERROR;
  mTCLResult.clear();
  if (mImpl->script) {
    // objv[0] is the name of this client
    if (COLVARSCRIPT_OK == mImpl->script->run(objc-1, objv+1)) {
      error_code = TCL_OK;
      mTCLResult = mImpl->get_error_msgs() + mImpl->script->str_result();
    }
  }
  // Cleanups
  for (int i = 0; i < objc; ++i) {
    delete[] objv[i];
  }
  delete[] objv;
  return error_code;
}

#endif // defined (NAMD_CUDA) || defined (NAMD_HIP)
