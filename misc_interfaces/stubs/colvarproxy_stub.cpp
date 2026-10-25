// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvars_version.h"

#define COLVARPROXY_VERSION COLVARS_VERSION


colvarproxy_stub::colvarproxy_stub()
{
  version_int = get_version_from_string(COLVARPROXY_VERSION);

  use_internal_pbc_ = true;

  b_simulation_running = false;

  // both fields are taken from data structures already available
  updated_masses_ = updated_charges_ = true;

  cvmodule = new colvarmodule(this);
  cvmodule->log("Using minimal testing interface.\n");

  cvmodule->cv_traj_freq = 0; // I/O will be handled explicitly
  cvmodule->restart_out_freq = 0;
  cvm::rotation::monitor_crossings = false; // Avoid unnecessary error messages

  cvmodule->setup_input();
  cvmodule->setup_output();

  colvarproxy_stub::setup();
}


colvarproxy_stub::~colvarproxy_stub()
{}


int colvarproxy_stub::setup()
{
  // Set system boundaries as the default (non-periodic)
  boundaries_.reset();

  if (cvmodule) {
    cvmodule->it = cvmodule->it_restart = 0;
    return cvmodule->update_engine_parameters();
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
  return true;
}


int colvarproxy_stub::set_unit_system(std::string const &units_in,
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


int colvarproxy_stub::read_frame_xyz(const char *filename, const bool write_force_file,
                                     int num_trial_evaluations,
                                     bool check_trial_consistency)
{
  int err = cvmodule->load_coords_xyz(filename, modify_atom_positions(), nullptr, true);
  std::string force_filename;
  if (write_force_file) {
    force_filename = cvmodule->output_prefix() + "_forces_" + cvm::to_str(cvmodule->it) + ".dat";
  }
  if ( !err ) {
    auto& new_forces = *(modify_atom_applied_forces());

    // The first frame must be a regular (committed) evaluation
    bool const trials_allowed = (cvmodule->step_relative() > 0);

    if (trials_allowed && (num_trial_evaluations > 0)) {
      std::vector<cvm::atom_pos> const saved_positions = *(modify_atom_positions());
      for (int k = 0; k < num_trial_evaluations; k++) {
        cvm::real const scale = 1.0 + 0.002 * cvm::real(k + 1);
        auto &positions = *(modify_atom_positions());
        for (size_t i = 0; i < positions.size(); i++) {
          positions[i] = saved_positions[i] * scale;
        }
        std::fill(new_forces.begin(), new_forces.end(), cvm::rvector{0, 0, 0});
        int const trial_err = cvmodule->calc_trial();
        if (trial_err != COLVARS_OK) return trial_err;
      }
      *(modify_atom_positions()) = saved_positions;
    }

    std::vector<cvm::rvector> trial_forces;
    if (trials_allowed && check_trial_consistency) {
      std::fill(new_forces.begin(), new_forces.end(), cvm::rvector{0, 0, 0});
      int const trial_err = cvmodule->calc_trial();
      if (trial_err != COLVARS_OK) return trial_err;
      trial_forces = new_forces;
    }

    std::fill(new_forces.begin(), new_forces.end(), cvm::rvector{0, 0, 0});
    cvmodule->calc();

    if (!trial_forces.empty()) {
      for (size_t i = 0; i < new_forces.size(); i++) {
        cvm::rvector const d = new_forces[i] - trial_forces[i];
        cvm::real const tol = 1.0e-10 * std::max(1.0, new_forces[i].norm());
        if (d.norm() > tol) {
          cvmodule->error("Error: trial and committed forces differ for atom " +
                          cvm::to_str(i) + ": " + cvm::to_str(trial_forces[i]) +
                          " vs. " + cvm::to_str(new_forces[i]) + "\n", COLVARS_BUG_ERROR);
          return COLVARS_BUG_ERROR;
        }
      }
    }

    cvmodule->it++;
    if (write_force_file) {
      const size_t numAtoms = modify_atom_positions()->size();
      std::ofstream ofs(force_filename);
      for (size_t i = 0; i < numAtoms; ++i) {
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << new_forces[i].x << std::setw(0) << " ";
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << new_forces[i].y << std::setw(0) << " ";
        ofs << std::scientific << std::setprecision(12) << std::setw(20) << new_forces[i].z << std::setw(0) << " ";
        ofs << std::endl;
      }
    }
  }
  return err;
}
