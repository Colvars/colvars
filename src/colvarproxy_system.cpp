// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <algorithm>
#include <cmath>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarproxy_system.h"

bool cvm::system_boundary_conditions::set_triclinic_minimum_image(cvm::rvector const &A,
                                                                 cvm::rvector const &B,
                                                                 cvm::rvector const &C)
{
  set_boundaries(true, true, true, A, B, C);
  if (type() == types::pbc_triclinic) {
    // The point initially lies in the centered fractional unit cell. Its
    // radius is bounded by the farthest corner. Any shorter image has each
    // fractional coordinate bounded by that radius times the dual norm.
    // A lattice shift outside the legacy 27 neighbors needs at least 1.5.
    const double radius2 = 0.25 * std::max({(A+B+C).norm2(), (A+B-C).norm2(),
                                           (A-B+C).norm2(), (A-B-C).norm2()});
    const double dual2 = std::max({reciprocal_cell_x.norm2(), reciprocal_cell_y.norm2(),
                                   reciprocal_cell_z.norm2()});
    if (radius2 * dual2 < 2.25) return true;
  }
  // Three-dimensional LLL reduction (delta=3/4). Integer column operations
  // preserve the lattice, including for long, strongly tilted cells.
  cvm::rvector basis[3] = {A, B, C}, orth[3], dual[3];
  int k = 1;
  for (int iteration = 0; iteration < 256; ++iteration) {
    for (int i = 0; i < 3; ++i) {
      orth[i] = basis[i];
      for (int j = 0; j < i; ++j)
        orth[i] -= (orth[i] * dual[j]) * orth[j];
      const double norm2 = orth[i].norm2();
      if (!(norm2 > 0.0) || !std::isfinite(norm2)) return false;
      dual[i] = orth[i] / norm2;
    }
    if (k == 3) {
      unit_cell_x = basis[0]; unit_cell_y = basis[1]; unit_cell_z = basis[2];
      reciprocal_cell_x = dual[0]; reciprocal_cell_y = dual[1]; reciprocal_cell_z = dual[2];
      periodic_x = periodic_y = periodic_z = true;
      set_type(types::pbc_reduced);
      return true;
    }
    for (int j = k-1; j >= 0; --j)
      basis[k] -= std::round(basis[k] * dual[j]) * basis[j];
    const double mu = basis[k] * dual[k-1];
    if (orth[k].norm2() >= (0.75 - mu*mu) * orth[k-1].norm2()) {
      ++k;
    } else {
      std::swap(basis[k], basis[k-1]);
      k = std::max(1, k-1);
    }
  }
  return false; // Numerically singular or insufficient precision for reduction.
}

colvarproxy_system::colvarproxy_system()
{
  angstrom_value_ = 0.0;
  kcal_mol_value_ = 0.0;
  timestep_ = 1.0;
  target_temperature_ = 0.0;
  boltzmann_ = 0.001987191; // Default: kcal/mol/K
  total_force_requested = false;
  indirect_lambda_biasing_force = 0.0;
  cached_alch_lambda_changed = false;
  cached_alch_lambda = -1.0;
}


colvarproxy_system::~colvarproxy_system() {}


int colvarproxy_system::set_unit_system(std::string const & /* units */,
                                        bool /* check_only */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_system::set_target_temperature(cvm::real T)
{
  target_temperature_ = T;
  return COLVARS_OK;
}


int colvarproxy_system::set_integration_timestep(cvm::real dt)
{
  timestep_ = dt;
  return COLVARS_OK;
}

int colvarproxy_system::set_time_step_factor(int fact)
{
  time_step_factor_ = fact;
  return COLVARS_OK;
}

cvm::real colvarproxy_system::rand_gaussian()
{
  // TODO define, document and implement a user method to set the value of this
  return 0.0;
}


void colvarproxy_system::add_energy(cvm::real /* energy */) {}


void colvarproxy_system::request_total_force(bool yesno)
{
  if (yesno == true)
    cvm::error_static("Error: total forces are currently not implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
}


bool colvarproxy_system::total_forces_enabled() const
{
  return false;
}


bool colvarproxy_system::total_forces_same_step() const
{
  return false;
}


cvm::rvector colvarproxy_system::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2) const
{
  return boundaries_.position_distance(pos1, pos2);
}


int colvarproxy_system::get_molid(int &)
{
  cvm::error_static("Error: only VMD allows the use of multiple \"molecules\", "
             "i.e. multiple molecular systems.", COLVARS_NOT_IMPLEMENTED);
  return -1;
}


int colvarproxy_system::get_alch_lambda(cvm::real * /* lambda */)
{
  return cvm::error_static("Error in get_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


void colvarproxy_system::set_alch_lambda(cvm::real lambda)
{
  cached_alch_lambda = lambda;
  cached_alch_lambda_changed = true;
}


int colvarproxy_system::send_alch_lambda()
{
  return cvm::error_static("Error in set_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::get_dE_dlambda(cvm::real * /* force */)
{
  return cvm::error_static("Error in get_dE_dlambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::apply_force_dE_dlambda(cvm::real* /* force */)
{
  return cvm::error_static("Error in apply_force_dE_dlambda: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::get_d2E_dlambda2(cvm::real*)
{
  return cvm::error_static("Error in get_d2E_dlambda2: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}
