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


cvm::rvector colvarproxy::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2) const
{
  return boundaries_.position_distance(pos1, pos2);
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
