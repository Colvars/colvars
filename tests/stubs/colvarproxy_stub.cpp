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
  if (colvars) {
    return colvars->setup();
  }
  return COLVARS_OK;
}


int colvarproxy_stub::set_unit_system(std::string const &units_in, 
                                            bool check_only)
{
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
      atoms_ncopies[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);

  if (aid < 0) {
    return INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);

  return index;
}

