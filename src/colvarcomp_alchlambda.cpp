// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include <algorithm>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"


colvar::alch_lambda::alch_lambda()
{
  set_function_type("alchLambda");

  provide(f_cvc_explicit_gradient, false);
  provide(f_cvc_gradient, false); // Cannot apply forces on this CVC
  provide(f_cvc_collect_atom_ids, false);

  provide(f_cvc_inv_gradient); // Projected force is TI derivative
  provide(f_cvc_Jacobian);     // Zero

  x.type(colvarvalue::type_scalar);

  // We need calculation every time step
  // default in Tinker-HP and NAMD2, must be enforced in NAMD3
  // Also checks back-end settings, ie. that alchemy is enabled
  // (in NAMD3: alchType TI, computeEnergies at the right frequency)
  cvm::proxy->request_alch_energy_freq(1);
  // TODO examine how this breaks everything - whereas alchOutFreq seems to work

  // Query initial value from back-end; will be overwritten if restarting from a state file
  cvm::proxy->get_alch_lambda(&x.real_value);
}


void colvar::alch_lambda::calc_value()
{
  // By default, follow external parameter
  // This might get overwritten by driving extended dynamics
  // (in apply_force() below)
  cvm::proxy->get_alch_lambda(&x.real_value);

  cvm::proxy->get_dE_dlambda(&ft.real_value);
  ft.real_value *= -1.0; // Convert energy derivative to force

  // Include any force due to bias on Flambda
  ft.real_value += cvm::proxy->indirect_lambda_biasing_force;
  cvm::proxy->indirect_lambda_biasing_force = 0.0;
}


void colvar::alch_lambda::calc_force_invgrads()
{
  // All the work is done in calc_value()
}


void colvar::alch_lambda::calc_Jacobian_derivative()
{
  jd = 0.0;
}


void colvar::alch_lambda::apply_force(colvarvalue const & /* force */)
{
  // new value will be cached and sent at end of timestep
  cvm::proxy->set_alch_lambda(x.real_value);
}


colvar::alch_Flambda::alch_Flambda()
{
  set_function_type("alch_Flambda");

  disable(f_cvc_explicit_gradient);
  disable(f_cvc_gradient);

  x.type(colvarvalue::type_scalar);
}


void colvar::alch_Flambda::calc_value()
{
  // Special workflow:
  // at the beginning of the timestep we get a force instead of calculating the value

  // Query initial value from back-end
  cvm::proxy->get_dE_dlambda(&x.real_value);
  x.real_value *= -1.0; // Energy derivative to force
}


void colvar::alch_Flambda::calc_gradients()
{
}


void colvar::alch_Flambda::apply_force(colvarvalue const &force)
{
  // Convert force on Flambda to force on dE/dlambda
  cvm::real f = -1.0 * force.real_value;
  // Send scalar force to back-end, which will distribute it onto atoms
  cvm::proxy->apply_force_dE_dlambda(&f);

  // Propagate force on Flambda to lambda internally
  cvm::real d2E_dlambda2;
  cvm::proxy->get_d2E_dlambda2(&d2E_dlambda2);

  // This accumulates a force, it needs to be zeroed when taken into account by lambda
  cvm::proxy->indirect_lambda_biasing_force += d2E_dlambda2 * f;
}

