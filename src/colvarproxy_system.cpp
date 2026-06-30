// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <algorithm>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarproxy.h"


int colvarproxy::set_unit_system(std::string const & /* units */,
                                        bool /* check_only */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy::set_target_temperature(cvm::real T)
{
  target_temperature_ = T;
  return COLVARS_OK;
}


int colvarproxy::set_integration_timestep(cvm::real dt)
{
  timestep_ = dt;
  return COLVARS_OK;
}

int colvarproxy::set_time_step_factor(int fact)
{
  time_step_factor_ = fact;
  return COLVARS_OK;
}

cvm::real colvarproxy::rand_gaussian()
{
  // TODO define, document and implement a user method to set the value of this
  return 0.0;
}


void colvarproxy::add_energy(cvm::real /* energy */) {}


void colvarproxy::request_total_force(bool yesno)
{
  if (yesno == true)
    cvm::error_static(cvmodule, "Error: total forces are currently not implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
}


bool colvarproxy::total_forces_enabled() const
{
  return false;
}


bool colvarproxy::total_forces_same_step() const
{
  return false;
}


inline int round_to_integer(cvm::real x)
{
  return int(cvm::floor(x+0.5));
}


void colvarproxy::update_pbc_lattice()
{
  // Periodicity is assumed in all directions

  if (boundaries_type == boundaries_unsupported ||
      boundaries_type == boundaries_non_periodic) {
    cvm::error_static(cvmodule, "Error: setting PBC lattice with unsupported boundaries.\n",
               COLVARS_BUG_ERROR);
    return;
  }

  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_y, unit_cell_z);
    reciprocal_cell_x = v/(v*unit_cell_x);
  }
  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_z, unit_cell_x);
    reciprocal_cell_y = v/(v*unit_cell_y);
  }
  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_x, unit_cell_y);
    reciprocal_cell_z = v/(v*unit_cell_z);
  }
}


void colvarproxy::reset_pbc_lattice()
{
  unit_cell_x.reset();
  unit_cell_y.reset();
  unit_cell_z.reset();
  reciprocal_cell_x.reset();
  reciprocal_cell_y.reset();
  reciprocal_cell_z.reset();
}


cvm::rvector colvarproxy::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2) const
{
  return position_distance_internal(pos1, pos2);
}


int colvarproxy::get_molid(int &)
{
  cvm::error_static(cvmodule, "Error: only VMD allows the use of multiple \"molecules\", "
             "i.e. multiple molecular systems.", COLVARS_NOT_IMPLEMENTED);
  return -1;
}


int colvarproxy::get_alch_lambda(cvm::real * /* lambda */)
{
  return cvm::error_static(cvmodule, "Error in get_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


void colvarproxy::set_alch_lambda(cvm::real lambda)
{
  cached_alch_lambda = lambda;
  cached_alch_lambda_changed = true;
}


int colvarproxy::send_alch_lambda()
{
  return cvm::error_static(cvmodule, "Error in set_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy::get_dE_dlambda(cvm::real * /* force */)
{
  return cvm::error_static(cvmodule, "Error in get_dE_dlambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy::apply_force_dE_dlambda(cvm::real* /* force */)
{
  return cvm::error_static(cvmodule, "Error in apply_force_dE_dlambda: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy::get_d2E_dlambda2(cvm::real*)
{
  return cvm::error_static(cvmodule, "Error in get_d2E_dlambda2: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}
