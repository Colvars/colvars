// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvaratoms.h"
#include "colvarproxy.h"

#include "colvarproxy_stub.h"


colvarproxy_stub::colvarproxy_stub()
{
  version_int = get_version_from_string(COLVARPROXY_VERSION);
  b_simulation_running = false;

  // both fields are taken from data structures already available
  updated_masses_ = updated_charges_ = true;

  colvars = new colvarmodule(this);
  cvm::log("Using minimal testing interface.\n");

  colvars->cv_traj_freq = 0; // I/O will be handled explicitly
  colvars->restart_out_freq = 0;
  cvm::rotation::monitor_crossings = false; // Avoid unnecessary error messages

  colvars->setup_input();
  colvars->setup_output();

  colvarproxy_stub::setup();
}


colvarproxy_stub::~colvarproxy_stub()
{}


int colvarproxy_stub::setup()
{
  boundaries_type = boundaries_non_periodic;
  reset_pbc_lattice();
  colvars->it = colvars->it_restart = 0;

  if (colvars) {
    return colvars->update_engine_parameters();
  }
  return COLVARS_OK;
}


void colvarproxy_stub::request_total_force(bool yesno)
{
  total_force_requested = yesno;
}


bool colvarproxy_stub::total_forces_enabled() const
{
  return total_force_requested;
}


bool colvarproxy_stub::total_forces_same_step() const
{
  return total_force_requested;
}


int colvarproxy_stub::set_unit_system(std::string const &units_in,
                                            bool check_only)
{
  // if check_only is specified, just test for compatibility
  // cvolvarmodule does that if new units are requested while colvars are already defined
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


void colvarproxy_stub::log(std::string const &message)
{
  std::cout << "colvars: " << message;
}


void colvarproxy_stub::error(std::string const &message)
{
  add_error_msg(message);
  std::cerr << "colvars: " << message;
}


int colvarproxy_stub::check_atom_id(int atom_number)
{
  return atom_number-1;
}


int colvarproxy_stub::init_atom(int atom_number)
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

  return index;
}


int colvarproxy_stub::read_frame_xyz(const char *filename)
{
  int err = colvars->load_coords_xyz(filename, modify_atom_positions(), nullptr, true);
  if ( !err ) {
    colvars->it++;
    colvars->calc();
  }
  return err;
}