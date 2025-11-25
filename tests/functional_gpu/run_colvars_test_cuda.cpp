// #define COLVARS_DEBUG true
#include "colvar_gpu_support.h"
#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvarproxy.h"
#include "CLI11.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#define COLVARPROXY_VERSION COLVARS_VERSION

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
    colvars = new colvarmodule(this);
    cvm::log("Using minimal CUDA testing interface.\n");

    colvars->cv_traj_freq = 0; // I/O will be handled explicitly
    colvars->restart_out_freq = 0;
    cvm::rotation::monitor_crossings = false; // Avoid unnecessary error messages

    colvars->setup_input();
    colvars->setup_output();

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
};

colvarproxy_stub_gpu::colvarproxy_stub_gpu():
  stream(0), d_mPositions(nullptr), d_mAppliedForces(nullptr),
  d_mTotalForces(nullptr), mAtomsChanged(false) {
  smp_mode = smp_mode_t::none;
  support_gpu = true;
  version_int = get_version_from_string(COLVARPROXY_VERSION);
  b_simulation_running = false;

  // both fields are taken from data structures already available
  updated_masses_ = updated_charges_ = true;

  checkGPUError(cudaStreamCreate(&stream));
}

colvarproxy_stub_gpu::~colvarproxy_stub_gpu() {
  checkGPUError(cudaStreamSynchronize(stream));
  checkGPUError(cudaStreamDestroy(stream));
  deallocateDeviceArrays();
}

int colvarproxy_stub_gpu::setup() {
  boundaries_type = boundaries_non_periodic;
  reset_pbc_lattice();
  colvars->it = colvars->it_restart = 0;
  if (colvars) {
    return colvars->update_engine_parameters();
  }
  return COLVARS_OK;
}

void colvarproxy_stub_gpu::request_total_force(bool yesno)
{
  total_force_requested = yesno;
}

bool colvarproxy_stub_gpu::total_forces_enabled() const
{
  return total_force_requested;
}

bool colvarproxy_stub_gpu::total_forces_same_step() const
{
  return total_force_requested;
}


int colvarproxy_stub_gpu::set_unit_system(std::string const &units_in,
                                            bool check_only)
{
  // if check_only is specified, just test for compatibility
  // colvarmodule sets this flag if new units are requested while colvars are already defined
  if (check_only) {
    if ((units != "" && units_in != units) || (units == "" && units_in != "real")) {
      cvm::error("Specified unit system \"" + units_in + "\" is incompatible with previous setting \""
                  + units + "\".\nReset the Colvars Module or delete all variables to change the unit.\n");
      return COLVARS_ERROR;
    } else {
      return COLVARS_OK;
    }
  }

  if (units_in == "real") {
    angstrom_value_ = 1.;
    kcal_mol_value_ = 1.;
  } else if (units_in == "metal") {
    angstrom_value_ = 1.;
    kcal_mol_value_ = 0.0433641017; // eV
    // inverse of LAMMPS value is 1/23.060549 = .043364102
  } else if (units_in == "electron") {
    angstrom_value_ = 1.88972612;    // Bohr
    kcal_mol_value_ = 0.00159360144; // Hartree
  } else if (units_in == "gromacs") {
    angstrom_value_ = 0.1;    // nm
    kcal_mol_value_ = 4.184;  // kJ/mol
  } else {
    cvm::error("Unknown unit system specified: \"" + units_in + "\". Supported are real, metal, electron, and gromacs.\n");
    return COLVARS_ERROR;
  }

  units = units_in;
  return COLVARS_OK;
}


void colvarproxy_stub_gpu::log(std::string const &message)
{
  std::cout << "colvars: " << message;
}


void colvarproxy_stub_gpu::error(std::string const &message)
{
  add_error_msg(message);
  std::cerr << "colvars: " << message;
  throw;
}


int colvarproxy_stub_gpu::check_atom_id(int atom_number)
{
  return atom_number-1;
}


int colvarproxy_stub_gpu::init_atom(int atom_number)
{
  // save time by checking first whether this atom has been requested before
  // (this is more common than a non-valid atom number)
  int aid = (atom_number-1);

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_refcount[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);

  if (aid < 0) {
    return COLVARS_INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);
  mAtomsChanged = true;
  return index;
}

void colvarproxy_stub_gpu::clear_atom(int index) {
  colvarproxy::clear_atom(index);
  mAtomsChanged = true;
}

int colvarproxy_stub_gpu::read_frame_xyz(const char *filename, const bool write_force_file)
{
  std::vector<cvm::rvector> positions(atoms_ids.size());
  int err = colvars->load_coords_xyz(filename, &positions, nullptr, true);
  // Convert to SOA and copy to GPU
  colvarproxy_atoms::atom_buffer_real_t positions_soa;
  const size_t numAtoms = positions.size();
  // if (numAtoms != positions.size()) {
  //   return cvm::error("Number of atoms mismatch!\n", COLVARS_ERROR);
  // }
  if (mAtomsChanged) {
    this->reallocate();
    if (colvars->gpu_calc) {
      // Need to rebuild the graph in case of reallocation
      colvars->gpu_calc->init();
    }
  }
  positions_soa.resize(3 * numAtoms);
  for (size_t i = 0; i < numAtoms; ++i) {
    positions_soa[i] = positions[i].x;
    positions_soa[i+numAtoms] = positions[i].y;
    positions_soa[i+2*numAtoms] = positions[i].z;
  }
  copy_HtoD(positions_soa.data(), d_mPositions, 3 * numAtoms);
  clear_device_array(d_mAppliedForces, 3 * numAtoms);
  clear_device_array(d_mTotalForces, 3 * numAtoms);
  std::string force_filename;
  if (write_force_file) {
    force_filename = colvars->output_prefix() + "_forces_" + cvm::to_str(colvars->it) + ".dat";
  }
  if ( !err ) {
    colvars->calc();
    colvars->it++;
    colvarproxy_atoms::atom_buffer_real_t h_applied_forces(3 * numAtoms);
    copy_DtoH(d_mAppliedForces, h_applied_forces.data(), 3 * numAtoms);
    if (write_force_file) {
      std::ofstream ofs(force_filename);
      for (size_t i = 0; i < numAtoms; ++i) {
        const cvm::real fx = h_applied_forces[i];
        const cvm::real fy = h_applied_forces[i + numAtoms];
        const cvm::real fz = h_applied_forces[i + 2 * numAtoms];
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << fx << std::setw(0) << " ";
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << fy << std::setw(0) << " ";
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << fz << std::setw(0) << " ";
        ofs << std::endl;
      }
    }
  }
  mAtomsChanged = false;
  return err;
}
#endif

int main(int argc, char *argv[]) {
#if defined (COLVARS_CUDA)
  CLI::App app{"Colvars stub interface for testing"};
  argv = app.ensure_utf8(argv);
  std::string configuration_file;
  std::string output_prefix;
  std::string trajectory_file;
  bool output_force = false;
  app.add_option("-c,--configuration_file", configuration_file, "Input Colvars configuration file")->required(true);
  app.add_option("-t,--trajectory_file", trajectory_file, "Input trajectory file")->required(true);
  const auto output_prefix_option =
    app.add_option("-o,--output_prefix", output_prefix, "Output file prefix")->required(false);
  app.add_flag("--force,!--no-force", output_force, "Write the force files")->needs(output_prefix_option);
  CLI11_PARSE(app, argc, argv);

  int err = COLVARS_OK;
  colvarproxy_stub_gpu *proxy = new colvarproxy_stub_gpu();
  proxy->init_cvm();
  // Initialize simple unit system to test file input
  err |= proxy->set_unit_system("real", false);
  // Initialize simple unit system to test file input
  err |= proxy->set_unit_system("real", false);

  if (argc > 3) {
    err |= proxy->set_output_prefix(output_prefix);
  }
  err |= proxy->colvars->setup_input();
  err |= proxy->colvars->setup_output();
  err |= proxy->colvars->read_config_file(configuration_file.c_str());
  if (err != COLVARS_OK) {
    cvm::log("Error occurred!\n");
  }

  if (argc > 2) {
    // Read number of atoms from XYZ header
    std::ifstream ifs(trajectory_file);
    int natoms;
    ifs >> natoms;
    ifs.close();
    cvm::log("Reading trajectory for " + cvm::to_str(natoms)
              + " atoms from XYZ file " + trajectory_file);
    for (int ai = 0; ai < natoms; ai++) {
      proxy->init_atom(ai+1);
    }
    err = cvm::get_error();
    if (err != COLVARS_OK) {
      cvm::log("Error occurred!\n");
    }
    int io_err = 0;
    while (!io_err) {
      io_err = proxy->read_frame_xyz(trajectory_file.c_str(), output_force);
      err = cvm::get_error();
      if (err != COLVARS_OK) {
        cvm::log("Error occurred!\n");
      }
      if (!io_err) cvm::log("Frame " + cvm::to_str(cvm::step_absolute()));
    }
    proxy->post_run();
    cvm::log("Done");
    err = cvm::get_error();
    if (err != COLVARS_OK) {
      cvm::log("Error occurred!\n");
    }
  }

  cvm::log("Input files read during this test:");
  unsigned char * args[2] = {
    (unsigned char *) "cv",
    (unsigned char *) "listinputfiles" };
  err |= run_colvarscript_command(2, args);
  cvm::log("  " + std::string(get_colvarscript_result()));

  double const max_gradient_error = proxy->colvars->get_max_gradient_error();
  if (max_gradient_error > 0.) {
    cvm::log("Max gradient error (debugGradients): " + cvm::to_str(max_gradient_error));

    double threshold = 1e-3;
    // Fail test if error is above threshold
    if (max_gradient_error > threshold) {
      cvm::log("Error: gradient inaccuracy is above threshold (" + cvm::to_str(threshold) + ")");
      err = 1;
    }
  }

  delete proxy;
  return err;
#else
  std::cout << "This program requires CUDA to test." << std::endl;
  return 1;
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
}
