#include "CudaGlobalMasterClient.h"
#include "colvarproxy_cudaglobalmaster.h"
#include "colvarproxy.h"
#include "Molecule.h"
#include "InfoStream.h"
#include "fstream_namd.h"
#include "CudaUtils.h"
#include "PDB.h"
#include "colvarparse.h"
#include "colvaratoms.h"

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

// Copied from colvarproxy_namd.C
enum e_pdb_field {
  e_pdb_none,
  e_pdb_occ,
  e_pdb_beta,
  e_pdb_x,
  e_pdb_y,
  e_pdb_z,
  e_pdb_ntot
};

// Copied from colvarproxy_namd.C
e_pdb_field pdb_field_str2enum(std::string const &pdb_field_str)
{
  e_pdb_field pdb_field = e_pdb_none;

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("O")) {
    pdb_field = e_pdb_occ;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("B")) {
    pdb_field = e_pdb_beta;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("X")) {
    pdb_field = e_pdb_x;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("Y")) {
    pdb_field = e_pdb_y;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("Z")) {
    pdb_field = e_pdb_z;
  }

  if (pdb_field == e_pdb_none) {
    cvm::error("Error: unsupported PDB field, \""+
               pdb_field_str+"\".\n", COLVARS_INPUT_ERROR);
  }

  return pdb_field;
}

class colvarproxy_impl: public colvarproxy {
public:
  colvarproxy_impl(const SimParameters* s, const Molecule* m);
  virtual ~colvarproxy_impl();
  void add_energy(cvm::real energy) override { mBiasEnergy += energy; }
  void log(std::string const &message) override;
  void error(std::string const &message) override;
  int init_atom(int atom_number) override;
  void clear_atom(int index) override;
  int check_atom_id(int atom_number) override;
  bool total_forces_enabled() const override { return total_force_requested; };
  bool total_forces_same_step() const override { return false; };
  int setup() override;
  int set_unit_system(std::string const &units_in, bool check_only) override;
  int check_replicas_enabled() override;
  int replica_index() override;
  int num_replicas() override;
  void replica_comm_barrier() override;
  int replica_comm_recv(char* msg_data, int buf_len, int src_rep) override;
  int replica_comm_send(char* msg_data, int msg_len, int dest_rep) override;
  std::ostream &output_stream(std::string const &output_name,
                              std::string const description) override;
  int flush_output_stream(std::string const &output_name) override;
  int flush_output_streams() override;
  int close_output_stream(std::string const &output_name) override;
  int close_output_streams() override;
  int backup_file(char const *filename) override;
  void initialize_from_cudagm(
    CudaGlobalMasterColvars* client,
    const std::vector<std::string>& arguments,
    const int deviceID, cudaStream_t stream);
  const bool atomsChanged() const {return mAtomsChanged;}
  int load_atoms_pdb(char const *filename,
                     cvm::atom_group &atoms,
                     std::string const &pdb_field,
                     double const pdb_field_value) override;
  int load_coords_pdb(char const *filename,
                      std::vector<cvm::atom_pos> &pos,
                      const std::vector<int> &indices,
                      std::string const &pdb_field,
                      double const pdb_field_value) override;
  void calculate();
  friend class CudaGlobalMasterColvars;
private:
  void allocateDeviceArrays();
  void deallocateDeviceArrays();
  void allocateHostArrays();
  void deallocateHostArrays();
  int update_target_temperature();
  double* d_mPositions;
  double* d_mAppliedForces;
  double* d_mTotalForces;
  double* d_mLattice;
  float*  d_mMass;
  float*  d_mCharges;
  double* h_mPositions;
  double* h_mAppliedForces;
  double* h_mTotalForces;
  double* h_mLattice;
  float*  h_mMass;
  float*  h_mCharges;
  double mBiasEnergy;
  cudaStream_t mStream;
  bool mAtomsChanged;
  bool first_timestep;
  cvm::step_number previous_NAMD_step;
  std::vector<std::string> mConfigFiles;
  const SimParameters* simParams;
  const Molecule* molecule;
  CudaGlobalMasterColvars* mClient;
  int m_device_id;
};

colvarproxy_impl::colvarproxy_impl(
  const SimParameters* s, const Molecule* m): colvarproxy(),
  d_mPositions(nullptr), d_mAppliedForces(nullptr),
  d_mTotalForces(nullptr), d_mLattice(nullptr),
  d_mMass(nullptr), d_mCharges(nullptr),
  h_mPositions(nullptr), h_mAppliedForces(nullptr),
  h_mTotalForces(nullptr), h_mLattice(nullptr),
  h_mMass(nullptr), h_mCharges(nullptr),
  mBiasEnergy(0), mAtomsChanged(false),
  first_timestep(true), previous_NAMD_step(0),
  simParams(s), molecule(m) {
}

colvarproxy_impl::~colvarproxy_impl() {
  deallocateDeviceArrays();
  deallocateHostArrays();
}

int colvarproxy_impl::setup() {
  engine_name_ = "NAMD_CUDAGLOBALMASTER";
  // if (colvars != nullptr) delete colvars;
  colvars = new colvarmodule(this);
  cvm::log("Using " + engine_name_ + " interface, version "+
           cvm::to_str(0)+".\n");
  colvars->cite_feature("NAMD engine");
  colvars->cite_feature("Colvars-NAMD interface");
  for (auto it = mConfigFiles.begin(); it != mConfigFiles.end(); ++it) {
    add_config("configfile", *it);
  }
  update_target_temperature();
  colvarproxy::parse_module_config();
  int error_code = colvarproxy::setup();
  if (colvars->size() == 0) {
    // Module is empty, nothing to do
    return COLVARS_OK;
  }
  if (simParams->wrapAll) {
    log("Warning: enabling wrapAll can lead to inconsistent results "
        "for Colvars calculations: please disable wrapAll, "
        "as is the default option in NAMD.\n");
  }
  set_integration_timestep(simParams->dt);
  log("updating target temperature (T = "+
      cvm::to_str(target_temperature())+" K).\n");
  error_code |= colvars->update_engine_parameters();
  error_code |= colvars->setup_input();
  error_code |= colvars->setup_output();
  if (simParams->firstTimestep != 0) {
    colvars->set_initial_step(static_cast<cvm::step_number>(simParams->firstTimestep));
  }
  return error_code;
}

void colvarproxy_impl::initialize_from_cudagm(
  CudaGlobalMasterColvars* client,
  const std::vector<std::string>& arguments,
  const int deviceID, cudaStream_t stream) {
  mClient = client;
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
  setup();
  colvarproxy_io::set_output_prefix(std::string(simParams->outputFilename));
  colvarproxy_io::set_restart_output_prefix(std::string(simParams->restartFilename));
  colvarproxy_io::set_default_restart_frequency(simParams->restartFrequency);
}

// Copied from colvarproxy_namd.C
int colvarproxy_impl::update_target_temperature()
{
  int error_code = COLVARS_OK;
  if (simParams->rescaleFreq > 0) {
    error_code |= set_target_temperature(simParams->rescaleTemp);
  } else if (simParams->reassignFreq > 0) {
    error_code |= set_target_temperature(simParams->reassignTemp);
  } else if (simParams->langevinOn) {
    error_code |= set_target_temperature(simParams->langevinTemp);
  } else if (simParams->tCoupleOn) {
    error_code |= set_target_temperature(simParams->tCoupleTemp);
  } else if (simParams->loweAndersenOn) {
    error_code |= set_target_temperature(simParams->loweAndersenTemp);
  } else if (simParams->stochRescaleOn) {
    error_code |= set_target_temperature(simParams->stochRescaleTemp);
  } else {
    error_code |= set_target_temperature(0.0);
  }
  return error_code;
}

int colvarproxy_impl::init_atom(int atom_number) {
  int aid = atom_number;
  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_refcount[i] += 1;
      return i;
    }
  }
  aid = check_atom_id(atom_number);
  if (aid < 0) {
    return aid;
  }
  int const index = colvarproxy::add_atom_slot(aid);
  // update_atom_properties(aid);
  mAtomsChanged = true;
  // TODO: This is an overkill!!!
  // I expect Colvars can either
  // (i) inform me the changing of atoms is completed
  // or (ii) initialize multiple atoms at once.
  deallocateDeviceArrays();
  allocateDeviceArrays();
  deallocateHostArrays();
  allocateHostArrays();
  return index;
}

void colvarproxy_impl::clear_atom(int index) {
  colvarproxy::clear_atom(index);
  mAtomsChanged = true;
  // TODO: This is an overkill!!!
  // I expect Colvars can either
  // (i) inform me the changing of atoms is completed
  // or (ii) initialize multiple atoms at once.
  deallocateDeviceArrays();
  allocateDeviceArrays();
  deallocateHostArrays();
  allocateHostArrays();
}

// Copied from colvarproxy_namd.C
int colvarproxy_impl::check_atom_id(int atom_number) {
  int const aid = atom_number;
  if (cvm::debug())
    log("Adding atom "+cvm::to_str(atom_number)+
        " for collective variables calculation.\n");
  if ((aid < 0) || (aid >= molecule->numAtoms)) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }
  return aid;
}

// Copied from colvarproxy_namd.C
int colvarproxy_impl::set_unit_system(std::string const &units_in, bool /*check_only*/) {
  cvm::log("units_in = " + units_in + "\n");
  if (units_in != "real") {
    cvm::error("Error: Specified unit system \"" + units_in + "\" is unsupported in NAMD. Supported units are \"real\" (A, kcal/mol).\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}

// Copied from colvarproxy_namd.C
void colvarproxy_impl::log(std::string const &message) {
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line))
    iout << "colvars: " << line << "\n";
  iout << endi;
}

// Copied from colvarproxy_namd.C
void colvarproxy_impl::error(std::string const &message)
{
  log(message);
  switch (cvm::get_error()) {
  case COLVARS_FILE_ERROR:
    errno = EIO; break;
  case COLVARS_NOT_IMPLEMENTED:
    errno = ENOSYS; break;
  case COLVARS_MEMORY_ERROR:
    errno = ENOMEM; break;
  }
  char const *msg = "Error in the collective variables module "
    "(see above for details)";
  if (errno) {
    NAMD_err(msg);
  } else {
    NAMD_die(msg);
  }
}

// Copied from colvarproxy_namd.C
int colvarproxy_impl::load_coords_pdb(char const *pdb_filename,
                                      std::vector<cvm::atom_pos> &pos,
                                      const std::vector<int> &indices,
                                      std::string const &pdb_field_str,
                                      double const pdb_field_value)
{
  if (pdb_field_str.size() == 0 && indices.size() == 0) {
    cvm::error("Bug alert: either PDB field should be defined or list of "
               "atom IDs should be available when loading atom coordinates!\n", COLVARS_BUG_ERROR);
  }

  e_pdb_field pdb_field_index;
  bool const use_pdb_field = (pdb_field_str.size() > 0);
  if (use_pdb_field) {
    pdb_field_index = pdb_field_str2enum(pdb_field_str);
  }

  // next index to be looked up in PDB file (if list is supplied)
  std::vector<int>::const_iterator current_index = indices.begin();

  PDB *pdb = new PDB(pdb_filename);
  size_t const pdb_natoms = pdb->num_atoms();

  if (pos.size() != pdb_natoms) {

    bool const pos_allocated = (pos.size() > 0);

    size_t ipos = 0, ipdb = 0;
    for ( ; ipdb < pdb_natoms; ipdb++) {

      if (use_pdb_field) {
        // PDB field mode: skip atoms with wrong value in PDB field
        double atom_pdb_field_value = 0.0;

        switch (pdb_field_index) {
        case e_pdb_occ:
          atom_pdb_field_value = (pdb->atom(ipdb))->occupancy();
          break;
        case e_pdb_beta:
          atom_pdb_field_value = (pdb->atom(ipdb))->temperaturefactor();
          break;
        case e_pdb_x:
          atom_pdb_field_value = (pdb->atom(ipdb))->xcoor();
          break;
        case e_pdb_y:
          atom_pdb_field_value = (pdb->atom(ipdb))->ycoor();
          break;
        case e_pdb_z:
          atom_pdb_field_value = (pdb->atom(ipdb))->zcoor();
          break;
        default:
          break;
        }

        if ( (pdb_field_value) &&
             (atom_pdb_field_value != pdb_field_value) ) {
          continue;
        } else if (atom_pdb_field_value == 0.0) {
          continue;
        }

      } else {
        // Atom ID mode: use predefined atom IDs from the atom group
        if (((int) ipdb) != *current_index) {
          // Skip atoms not in the list
          continue;
        } else {
          current_index++;
        }
      }

      if (!pos_allocated) {
        pos.push_back(cvm::atom_pos(0.0, 0.0, 0.0));
      } else if (ipos >= pos.size()) {
        cvm::error("Error: the PDB file \""+
                   std::string(pdb_filename)+
                   "\" contains coordinates for "
                   "more atoms than needed.\n", COLVARS_BUG_ERROR);
      }

      pos[ipos] = cvm::atom_pos((pdb->atom(ipdb))->xcoor(),
                                (pdb->atom(ipdb))->ycoor(),
                                (pdb->atom(ipdb))->zcoor());
      ipos++;
      if (!use_pdb_field && current_index == indices.end())
        break;
    }

    if (ipos < pos.size() || (!use_pdb_field && current_index != indices.end())) {
      size_t n_requested = use_pdb_field ? pos.size() : indices.size();
      cvm::error("Error: number of matching records in the PDB file \""+
                 std::string(pdb_filename)+"\" ("+cvm::to_str(ipos)+
                 ") does not match the number of requested coordinates ("+
                 cvm::to_str(n_requested)+").\n", COLVARS_INPUT_ERROR);
      return COLVARS_ERROR;
    }
  } else {

    // when the PDB contains exactly the number of atoms of the array,
    // ignore the fields and just read coordinates
    for (size_t ia = 0; ia < pos.size(); ia++) {
      pos[ia] = cvm::atom_pos((pdb->atom(ia))->xcoor(),
                              (pdb->atom(ia))->ycoor(),
                              (pdb->atom(ia))->zcoor());
    }
  }

  delete pdb;
  return COLVARS_OK;
}

// Copied from colvarproxy_namd.C
int colvarproxy_impl::load_atoms_pdb(char const *pdb_filename,
                                     cvm::atom_group &atoms,
                                     std::string const &pdb_field_str,
                                     double const pdb_field_value)
{
  if (pdb_field_str.size() == 0)
    cvm::error("Error: must define which PDB field to use "
               "in order to define atoms from a PDB file.\n", COLVARS_INPUT_ERROR);

  PDB *pdb = new PDB(pdb_filename);
  size_t const pdb_natoms = pdb->num_atoms();

  e_pdb_field pdb_field_index = pdb_field_str2enum(pdb_field_str);

  for (size_t ipdb = 0; ipdb < pdb_natoms; ipdb++) {

    double atom_pdb_field_value = 0.0;

    switch (pdb_field_index) {
    case e_pdb_occ:
      atom_pdb_field_value = (pdb->atom(ipdb))->occupancy();
      break;
    case e_pdb_beta:
      atom_pdb_field_value = (pdb->atom(ipdb))->temperaturefactor();
      break;
    case e_pdb_x:
      atom_pdb_field_value = (pdb->atom(ipdb))->xcoor();
      break;
    case e_pdb_y:
      atom_pdb_field_value = (pdb->atom(ipdb))->ycoor();
      break;
    case e_pdb_z:
      atom_pdb_field_value = (pdb->atom(ipdb))->zcoor();
      break;
    default:
      break;
    }

    if ( (pdb_field_value) &&
         (atom_pdb_field_value != pdb_field_value) ) {
      continue;
    } else if (atom_pdb_field_value == 0.0) {
      continue;
    }

    if (atoms.is_enabled(colvardeps::f_ag_scalable)) {
      atoms.add_atom_id(ipdb);
    } else {
      atoms.add_atom(cvm::atom(ipdb+1));
    }
  }

  delete pdb;
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
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

void colvarproxy_impl::allocateHostArrays() {
  const int numAtoms = atoms_ids.size();
  allocate_host<double>(&h_mPositions, 3*numAtoms);
  allocate_host<double>(&h_mAppliedForces, 3*numAtoms);
  if (mClient->requestedTotalForcesAtomsChanged()) {
    allocate_host<double>(&h_mTotalForces, 3*numAtoms);
  }
  allocate_host<float>(&h_mMass, numAtoms);
  allocate_host<float>(&h_mCharges, numAtoms);
}

void colvarproxy_impl::deallocateHostArrays() {
  deallocate_host<double>(&h_mPositions);
  deallocate_host<double>(&h_mAppliedForces);
  if (mClient->requestedTotalForcesAtomsChanged()) {
    deallocate_host<double>(&h_mTotalForces);
  }
  deallocate_host<float>(&h_mMass);
  deallocate_host<float>(&h_mCharges);
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

void colvarproxy_impl::calculate() {
  const int64_t step = mClient->getStep();
  if (first_timestep) {
    // TODO: Do I really need to call them again?
    // setup();
    update_target_temperature();
    colvars->update_engine_parameters();
    colvars->setup_input();
    colvars->setup_output();
    first_timestep = false;
  } else {
    if ( step - previous_NAMD_step == 1 ) {
      colvars->it++;
      b_simulation_continuing = false;
    } else {
      // Cases covered by this condition:
      // - run 0
      // - beginning of a new run statement
      // The internal counter is not incremented, and the objects are made
      // aware of this via the following flag
      b_simulation_continuing = true;

      // Update NAMD output and restart prefixes
      colvarproxy_io::set_output_prefix(std::string(simParams->outputFilename));
      colvarproxy_io::set_restart_output_prefix(std::string(simParams->restartFilename));
      colvarproxy_io::set_default_restart_frequency(simParams->restartFrequency);
      colvars->setup_output();
    }
  }
  previous_NAMD_step = step;
  // Clear the previous bias energy
  mBiasEnergy = 0;
  // If the atom selection is changed, re-allocate arrays
  int savedDevice;
  cudaCheck(cudaGetDevice(&savedDevice));
  cudaCheck(cudaSetDevice(m_device_id));
  // TODO: Colvars does not support GPU, so we have to copy the buffers manually
  const size_t numAtoms = atoms_ids.size();
  copy_DtoH(d_mPositions, h_mPositions, 3*numAtoms, mStream);
  if (mClient->requestedTotalForcesAtomsChanged()) {
    copy_DtoH(d_mTotalForces, h_mTotalForces, 3*numAtoms, mStream);
  }
  if (mClient->requestUpdateMasses()) {
    copy_DtoH(d_mMass, h_mMass, numAtoms, mStream);
  }
  if (mClient->requestUpdateCharges()) {
    copy_DtoH(d_mCharges, h_mCharges, numAtoms, mStream);
  }
  if (mClient->requestUpdateLattice()) {
    copy_DtoH(d_mLattice, h_mLattice, 3*4, mStream);
  }
  // Synchronize the stream to make sure the host buffers are ready
  cudaCheck(cudaStreamSynchronize(mStream));
  // Transform the arrays for Colvars
  std::vector<cvm::atom_pos> &colvars_pos = *(modify_atom_positions());
  for (size_t i = 0; i < numAtoms; ++i) {
    colvars_pos[i] = {h_mPositions[i],
                      h_mPositions[i+numAtoms],
                      h_mPositions[i+2*numAtoms]};
  }
  if (mClient->requestedTotalForcesAtomsChanged()) {
    std::vector<cvm::rvector> &colvars_total_force = *(modify_atom_total_forces());
    for (size_t i = 0; i < numAtoms; ++i) {
      colvars_total_force[i] = {
        h_mTotalForces[i],
        h_mTotalForces[i+numAtoms],
        h_mTotalForces[i+2*numAtoms]};
    }
  }
  if (mClient->requestUpdateMasses()) {
    std::vector<cvm::real> &colvars_mass = *(modify_atom_masses());
    for (size_t i = 0; i < numAtoms; ++i) {
      colvars_mass[i] = cvm::real(h_mMass[i]);
    }
  }
  if (mClient->requestUpdateCharges()) {
    std::vector<cvm::real> &colvars_charge  = *(modify_atom_charges());
    for (size_t i = 0; i < numAtoms; ++i) {
      colvars_charge[i] = cvm::real(h_mCharges[i]);
    }
  }
  if (mClient->requestUpdateLattice()) {
    unit_cell_x.set(h_mLattice[0], h_mLattice[1], h_mLattice[2]);
    unit_cell_y.set(h_mLattice[3], h_mLattice[4], h_mLattice[5]);
    unit_cell_z.set(h_mLattice[6], h_mLattice[7], h_mLattice[8]);
    const Vector a1(h_mLattice[0], h_mLattice[1], h_mLattice[2]);
    const Vector a2(h_mLattice[3], h_mLattice[4], h_mLattice[5]);
    const Vector a3(h_mLattice[6], h_mLattice[7], h_mLattice[8]);
    const int p1 = ( a1.length2() ? 1 : 0 );
    const int p2 = ( a2.length2() ? 1 : 0 );
    const int p3 = ( a3.length2() ? 1 : 0 );
    if (!p1 && !p2 && !p3) {
      boundaries_type = boundaries_non_periodic;
      reset_pbc_lattice();
    } else if (p1 && p2 && p3) {
      if (( ! ( a1.y || a1.z || a2.x || a2.z || a3.x || a3.y ) )) {
        boundaries_type = boundaries_pbc_ortho;
      } else {
        boundaries_type = boundaries_pbc_triclinic;
      }
      colvarproxy_system::update_pbc_lattice();
    } else {
      boundaries_type = boundaries_unsupported;
    }
  }
  // Run Colvars
  if (cvm::debug()) {
    print_input_atomic_data();
  }
  if (colvars->calc() != COLVARS_OK) {
    cvm::error("Error in the collective variables module.\n", COLVARS_ERROR);
  }
  if (cvm::debug()) {
    print_output_atomic_data();
  }
  // Update applied forces
  std::vector<cvm::rvector> &colvars_applied_force
    = *(modify_atom_applied_forces());
  for (size_t i = 0; i < numAtoms; ++i) {
    h_mAppliedForces[i] = colvars_applied_force[i].x;
    h_mAppliedForces[i+numAtoms] = colvars_applied_force[i].y;
    h_mAppliedForces[i+2*numAtoms] = colvars_applied_force[i].z;
  }
  // TODO: Why do I need to clean the applied forces manually?
  std::fill(colvars_applied_force.begin(),
            colvars_applied_force.end(),
            cvm::rvector(0, 0, 0));
  // log("NUMBER OF ATOMS IN COLVARS: " + cvm::to_str(numAtoms) + "\n");
  copy_HtoD(h_mAppliedForces, d_mAppliedForces, 3*numAtoms, mStream);
  // NOTE: I think I can skip the syncrhonization here because this client
  //       share the same stream as the CudaGlobalMasterServer object
  // cudaCheck(cudaStreamSynchronize(mStream));
  // NAMD does not destruct GlobalMaster objects, so we must remember
  // to write all output files at the end of a run
  if (step == simParams->N) {
    post_run();
  }
  // Restore the GPU device
  cudaCheck(cudaSetDevice(savedDevice));
  // Update the mAtomsChanged flag
  mAtomsChanged = false;
}

std::ostream & colvarproxy_impl::output_stream(std::string const &output_name,
                                               std::string const description)
{
  if (cvm::debug()) {
    cvm::log("Using colvarproxy_namd::output_stream()\n");
  }

  if (!io_available()) {
    cvm::error("Error: trying to access an output file/channel "
               "from the wrong thread.\n", COLVARS_BUG_ERROR);
    return *output_stream_error_;
  }

  if (output_streams_.count(output_name) > 0) {
    return *(output_streams_[output_name]);
  }

  backup_file(output_name.c_str());

  output_streams_[output_name] = new ofstream_namd(output_name.c_str(), std::ios::binary);
  if (! output_streams_[output_name]->good()) {
    cvm::error("Error: cannot write to "+description+" \""+output_name+"\".\n",
               COLVARS_FILE_ERROR);
  }

  return *(output_streams_[output_name]);
}


int colvarproxy_impl::flush_output_stream(std::string const &output_name)
{
  if (!io_available()) {
    return COLVARS_OK;
  }

  if (output_streams_.count(output_name) > 0) {
    (reinterpret_cast<ofstream_namd *>(output_streams_[output_name]))->flush();
    return COLVARS_OK;
  }

  return cvm::error("Error: trying to flush an output file/channel "
                    "that wasn't open.\n", COLVARS_BUG_ERROR);
}


int colvarproxy_impl::flush_output_streams()
{
  if (!io_available()) {
    return COLVARS_OK;
  }

  for (std::map<std::string, std::ostream *>::iterator osi = output_streams_.begin();
       osi != output_streams_.end();
       osi++) {
    (reinterpret_cast<ofstream_namd *>(osi->second))->flush();
  }

  return COLVARS_OK;
}


int colvarproxy_impl::close_output_stream(std::string const &output_name)
{
  if (!io_available()) {
    return cvm::error("Error: trying to access an output file/channel "
                      "from the wrong thread.\n", COLVARS_BUG_ERROR);
  }

  if (output_streams_.count(output_name) > 0) {
    (reinterpret_cast<ofstream_namd *>(output_streams_[output_name]))->close();
    delete output_streams_[output_name];
    output_streams_.erase(output_name);
  }

  return COLVARS_OK;
}


int colvarproxy_impl::close_output_streams()
{
  if (! io_available()) {
    return COLVARS_OK;
  }

  for (std::map<std::string, std::ostream *>::iterator osi = output_streams_.begin();
       osi != output_streams_.end();
       osi++) {
    (reinterpret_cast<ofstream_namd *>(osi->second))->close();
  }
  output_streams_.clear();

  return COLVARS_OK;
}


int colvarproxy_impl::backup_file(char const *filename)
{
  if (std::string(filename).rfind(std::string(".colvars.state")) != std::string::npos) {
    NAMD_backup_file(filename, ".old");
  } else {
    NAMD_backup_file(filename, ".BAK");
  }
  return COLVARS_OK;
}

CudaGlobalMasterColvars::CudaGlobalMasterColvars():
  CudaGlobalMasterClient()
{
  mImpl = std::make_unique<colvarproxy_impl>(
    CudaGlobalMasterClient::getSimParameters(),
    CudaGlobalMasterClient::getMolecule());
}

CudaGlobalMasterColvars::~CudaGlobalMasterColvars() {}

void CudaGlobalMasterColvars::initialize(
  const std::vector<std::string>& arguments,
  int deviceID, cudaStream_t stream) {
  CudaGlobalMasterClient::initialize(arguments, deviceID, stream);
  mImpl->initialize_from_cudagm(this, arguments, deviceID, stream);
  // iout << "Using stream: " << stream << "\n" << endi;
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
