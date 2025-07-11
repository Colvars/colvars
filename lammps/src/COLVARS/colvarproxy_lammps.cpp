// clang-format off
// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarproxy_lammps.h"

#include "domain.h"
#include "error.h"
#include "force.h"
#include "lammps.h"             // includes <cstdio>, <mpi.h>, <string>, <vector>
#include "update.h"
#include "utils.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"

#define HASH_FAIL  -1

/* ---------------------------------------------------------------------- */

colvarproxy_lammps::colvarproxy_lammps(LAMMPS_NS::LAMMPS *lmp)  : _lmp(lmp), _random(nullptr)
{
  engine_name_ = "LAMMPS";

  first_timestep = true;
  previous_step = -1;
  do_exit = false;

  bias_energy = 0.0;

  engine_ready_ = false;
}

/* ---------------------------------------------------------------------- */

void colvarproxy_lammps::init()
{
  version_int = get_version_from_string(COLVARPROXY_VERSION);

  // create the colvarmodule instance
  colvars = new colvarmodule(this);

  // Create instance of scripting interface
  script = new colvarscript(this, colvars);

  cvm::log("Using LAMMPS interface, version " + cvm::to_str(COLVARPROXY_VERSION) + ".\n");

  colvars->cite_feature("LAMMPS engine");
  colvars->cite_feature("Colvars-LAMMPS interface");

  angstrom_value_ = _lmp->force->angstrom;
  boltzmann_ = _lmp->force->boltz;
  set_integration_timestep(_lmp->update->dt * _lmp->force->femtosecond);

  if (_lmp->update->ntimestep != 0) {
    colvars->set_initial_step(static_cast<cvm::step_number>(_lmp->update->ntimestep));
  }
}

/* ---------------------------------------------------------------------- */

colvarproxy_lammps::~colvarproxy_lammps()
{
  if (_random) delete _random;
}

/* ---------------------------------------------------------------------- */

void colvarproxy_lammps::set_random_seed(int seed)
{
  if (_random) delete _random;

  _random = new LAMMPS_NS::RanPark(_lmp, seed);
}

/* ----------------------------------------------------------------------
   re-initialize data where needed
------------------------------------------------------------------------- */

int colvarproxy_lammps::setup()
{
  int error_code = colvarproxy::setup();
  set_integration_timestep(_lmp->update->dt * _lmp->force->femtosecond);
  error_code |= colvars->update_engine_parameters();
  error_code |= colvars->setup_input();
  error_code |= colvars->setup_output();
  return error_code;
}

/* ----------------------------------------------------------------------
   trigger colvars computation
------------------------------------------------------------------------- */

double colvarproxy_lammps::compute()
{
  if (cvm::debug()) {
    cvm::log(std::string(cvm::line_marker) +
        "colvarproxy_lammps step no. " +
        cvm::to_str(_lmp->update->ntimestep) + " [first - last = " +
        cvm::to_str(_lmp->update->beginstep) + " - " +
        cvm::to_str(_lmp->update->endstep) + "]\n");
  }

  if (first_timestep) {
    first_timestep = false;
  } else {
    // Use the time step number from LAMMPS Update object
    if (_lmp->update->ntimestep - previous_step == 1) {
      colvarmodule::it++;
      b_simulation_continuing = false;
    } else {
      // Cases covered by this condition:
      // - run 0
      // - beginning of a new run statement
      // The internal counter is not incremented, and the objects are made
      // aware of this via the following flag
      b_simulation_continuing = true;
    }

  }
  previous_step = _lmp->update->ntimestep;

  unit_cell_x.set(_lmp->domain->xprd, 0.0, 0.0);
  unit_cell_y.set(0.0, _lmp->domain->yprd, 0.0);
  unit_cell_z.set(0.0, 0.0, _lmp->domain->zprd);

  if (_lmp->domain->xperiodic == 0 && _lmp->domain->yperiodic == 0 &&
      _lmp->domain->zperiodic == 0) {
    boundaries_type = boundaries_non_periodic;
    reset_pbc_lattice();
  } else if ((_lmp->domain->nonperiodic == 0) &&
             (_lmp->domain->dimension == 3) &&
             (_lmp->domain->triclinic == 0)) {
    // Orthogonal unit cell
    boundaries_type = boundaries_pbc_ortho;
    colvarproxy_system::update_pbc_lattice();
    // It is safer to let LAMMPS deal with high-tilt triclinic boxes
  } else {
    boundaries_type = boundaries_unsupported;
  }

  if (cvm::debug()) {
    cvm::log(std::string(cvm::line_marker) +
             "colvarproxy_lammps, step no. " + cvm::to_str(colvarmodule::it) + "\n" +
             "Updating internal data.\n");
  }

  // zero the forces on the atoms, so that they can be accumulated by the colvars
  for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++)
    atoms_new_colvar_forces[i].reset();

  bias_energy = 0.0;

  if (cvm::debug()) {
    cvm::log("atoms_ids = " + cvm::to_str(atoms_ids) + "\n");
    cvm::log("atoms_refcount = " + cvm::to_str(atoms_refcount) + "\n");
    cvm::log("atoms_positions = " + cvm::to_str(atoms_positions) + "\n");
    cvm::log("atoms_new_colvar_forces = " + cvm::to_str(atoms_new_colvar_forces) + "\n");
  }

  // Call the collective variable module
  if (colvars->calc() != COLVARS_OK)
    cvm::error("Error in the collective variables module.\n", COLVARS_ERROR);

  if (cvm::debug()) {
    cvm::log("atoms_ids = " + cvm::to_str(atoms_ids) + "\n");
    cvm::log("atoms_refcount = " + cvm::to_str(atoms_refcount) + "\n");
    cvm::log("atoms_positions = " + cvm::to_str(atoms_positions) + "\n");
    cvm::log("atoms_new_colvar_forces = " + cvm::to_str(atoms_new_colvar_forces) + "\n");
  }

  return bias_energy;
}

/* ---------------------------------------------------------------------- */

cvm::rvector colvarproxy_lammps::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2)
  const
{
  double xtmp = pos2.x - pos1.x;
  double ytmp = pos2.y - pos1.y;
  double ztmp = pos2.z - pos1.z;
  _lmp->domain->minimum_image_big(FLERR, xtmp,ytmp,ztmp);
  return {xtmp, ytmp, ztmp};
}

/* ---------------------------------------------------------------------- */

void colvarproxy_lammps::log(std::string const &message)
{
  LAMMPS_NS::utils::logmesg(_lmp, message);
}

/* ---------------------------------------------------------------------- */

void colvarproxy_lammps::error(std::string const &message)
{
  log(message);
  _lmp->error->one(FLERR, "Fatal error in the collective variables module");
}

/* ---------------------------------------------------------------------- */

char const *colvarproxy_lammps::script_obj_to_str(unsigned char *obj)
{
  // For now we assume that all objects passed by FixColvars are strings
  return reinterpret_cast<char *>(obj);
}

/* ---------------------------------------------------------------------- */

std::vector<std::string> colvarproxy_lammps::script_obj_to_str_vector(unsigned char *obj)
{
  if (cvm::debug()) {
    cvm::log("Called colvarproxy_lammps::script_obj_to_str_vector().\n");
  }
  std::string const input(reinterpret_cast<char *>(obj));
  return LAMMPS_NS::utils::split_words(input); // :-)))
}

/* ---------------------------------------------------------------------- */

int colvarproxy_lammps::set_unit_system(std::string const &units_in, bool /*check_only*/)
{
  std::string lmp_units = _lmp->update->unit_style;
  if (units_in != lmp_units) {
    cvm::error("Error: Specified unit system for Colvars \"" + units_in  +
               "\" is incompatible with LAMMPS internal units (" + lmp_units + ").\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}



int colvarproxy_lammps::check_atom_id(int atom_number)
{
  int const aid = atom_number;

  if (cvm::debug())
    log("Adding atom " + cvm::to_str(atom_number) + " for collective variables calculation.\n");

  // TODO add upper boundary check?
  if ((aid < 0)) {
    cvm::error("Error: invalid atom number specified, "  +
               cvm::to_str(atom_number) + "\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  return aid;
}

/* ---------------------------------------------------------------------- */

int colvarproxy_lammps::init_atom(int atom_number)
{
  int aid = atom_number;

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_refcount[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);
  if (aid < 0) return aid;

  int const index = colvarproxy::add_atom_slot(aid);
  // add entries for the LAMMPS-specific fields
  atoms_types.push_back(0);

  return index;
}
