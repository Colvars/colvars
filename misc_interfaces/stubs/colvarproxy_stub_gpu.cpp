// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include <iostream>
#include <fstream>
// #include <string>
#include <iomanip>

#include "colvar_gpu_support.h"
#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvarproxy.h"
#include "colvarproxy_stub_gpu.h"

#define COLVARPROXY_VERSION COLVARS_VERSION

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)

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
  int error_code = colvarproxy::setup();
  // Set system boundaries as the default (non-periodic)
  boundaries_.reset();
  if (cvmodule) {
    cvmodule->it = cvmodule->it_restart = 0;
    return cvmodule->update_engine_parameters();
  }
  return error_code;
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
      cvmodule->error("Specified unit system \"" + units_in + "\" is incompatible with previous setting \""
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
    cvmodule->error("Unknown unit system specified: \"" + units_in + "\". Supported are real, metal, electron, and gromacs.\n");
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
  // std::vector<cvm::rvector> positions(atoms_ids.size());
  m_positions.resize(atoms_ids.size());
  int error_code = cvmodule->load_coords_xyz(filename, &m_positions, nullptr, true);
  // Convert to SOA and copy to GPU
  const size_t numAtoms = m_positions.size();
  // cvmodule->log("numAtoms = " + cvm::to_str(numAtoms) + "\n");
  if (mAtomsChanged) {
    this->reallocate();
    error_code |= cvmodule->proxy_gpu_buffers_reallocated_done();
  }
  // The XYZ file lists atoms by increasing ID, but atom slots follow the order
  // in which the atoms were requested
  m_positions_soa.resize(3 * numAtoms);
  for (size_t i = 0; i < numAtoms; ++i) {
    m_positions_soa[i] = m_positions[atoms_ids[i]].x;
    m_positions_soa[i+numAtoms] = m_positions[atoms_ids[i]].y;
    m_positions_soa[i+2*numAtoms] = m_positions[atoms_ids[i]].z;
  }
  error_code |= copy_HtoD_async(m_positions_soa.data(), d_mPositions, 3 * numAtoms, stream);
  error_code |= clear_device_array_async(d_mAppliedForces, 3 * numAtoms, stream);
  error_code |= clear_device_array_async(d_mTotalForces, 3 * numAtoms, stream);
  error_code |= checkGPUError(cudaEventRecord(get_event(colvarproxy_gpu::event_type::copy_atoms), stream));
  std::string force_filename;
  if (write_force_file) {
    force_filename = cvmodule->output_prefix() + "_forces_" + cvm::to_str(cvmodule->it) + ".dat";
  }
  if ( !error_code ) {
    error_code |= cvmodule->calc();
    cvmodule->it++;
    m_applied_forces_soa.resize(3 * numAtoms);
    error_code |= copy_DtoH_async(d_mAppliedForces, m_applied_forces_soa.data(), 3 * numAtoms, stream);
    error_code |= checkGPUError(cudaStreamSynchronize(stream));
    if (write_force_file) {
      // Write the forces by increasing atom ID
      m_forces_by_id.resize(numAtoms);
      for (size_t i = 0; i < numAtoms; ++i) {
        m_forces_by_id[atoms_ids[i]] = cvm::rvector(m_applied_forces_soa[i],
                                                    m_applied_forces_soa[i + numAtoms],
                                                    m_applied_forces_soa[i + 2 * numAtoms]);
      }
      std::ofstream ofs(force_filename);
      for (size_t i = 0; i < numAtoms; ++i) {
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << m_forces_by_id[i].x << std::setw(0) << " ";
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << m_forces_by_id[i].y << std::setw(0) << " ";
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << m_forces_by_id[i].z << std::setw(0) << " ";
        ofs << std::endl;
      }
    }
  }
  mAtomsChanged = false;
  return error_code;
}
#endif

