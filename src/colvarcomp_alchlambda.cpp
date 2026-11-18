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
  init_scalar_boundaries(0.0, 1.0);

  // Query initial value from back-end; will be overwritten if restarting from a state file
  cvmodule->proxy->get_alch_lambda(&x.real_value);
}


int colvar::alch_lambda::init_alchemy(int factor)
{
  // We need calculation every time step
  // default in Tinker-HP and NAMD2, must be enforced in NAMD3
  // Also checks back-end settings, ie. that alchemy is enabled
  // (in NAMD3: alchType TI, computeEnergies at the right frequency)

  // Forbid MTS until fully implemented
  if (factor != 1) {
    return cvmodule->error("Error: timeStepFactor > 1 is not yet supported for alchemical variables.");
  }
  cvmodule->proxy->request_alch_energy_freq(factor);

  return COLVARS_OK;
}


void colvar::alch_lambda::calc_value()
{
  // By default, follow external parameter
  // This might get overwritten by driving extended dynamics
  // (in apply_force() below)
  cvmodule->proxy->get_alch_lambda(&x.real_value);

  cvmodule->proxy->get_dE_dlambda(&ft.real_value);
  ft.real_value *= -1.0; // Convert energy derivative to force

  // Include any force due to bias on Flambda
  ft.real_value += cvmodule->proxy->indirect_lambda_biasing_force;
  cvmodule->proxy->indirect_lambda_biasing_force = 0.0;
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
  // Forces, if any, are applied in colvar::update_extended_Lagrangian()
}


int colvar::alch_lambda::set_external_value(bool now)
{
  // Cache value to be communicated to back-end between time steps
  cvmodule->proxy->set_alch_lambda(x.real_value);
  if (now) {
    // If requested (e.g. upon restarting), sync to back-end
    return cvmodule->proxy->send_alch_lambda();
  }
  return COLVARS_OK;
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
  cvmodule->proxy->get_dE_dlambda(&x.real_value);
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
  cvmodule->proxy->apply_force_dE_dlambda(&f);

  // Propagate force on Flambda to lambda internally
  cvm::real d2E_dlambda2;
  cvmodule->proxy->get_d2E_dlambda2(&d2E_dlambda2);

  // This accumulates a force, it needs to be zeroed when taken into account by lambda
  cvmodule->proxy->indirect_lambda_biasing_force += d2E_dlambda2 * f;
}


colvar::tabf_alpha::tabf_alpha()
{
  set_function_type("tabfAlpha");

  provide(f_cvc_explicit_gradient, false);
  provide(f_cvc_gradient, false); // Cannot apply forces on this CVC
  provide(f_cvc_collect_atom_ids, false);

  provide(f_cvc_inv_gradient); // Projected force is the alpha derivative
  provide(f_cvc_Jacobian);     // Zero

  x.type(colvarvalue::type_scalar);
}


int colvar::tabf_alpha::init(std::string const &conf)
{
  int error_code = cvc::init(conf);
  if (error_code != COLVARS_OK) return error_code;

  // Use the back-end value until a restart supplies one.
  auto *proxy = cvmodule->proxy;
  error_code = proxy->get_tabf_alpha(&x.real_value);
  if (error_code != COLVARS_OK) return proxy->report_tabf_error(error_code, "get_tabf_alpha");

  cvm::real lower = 0.0;
  cvm::real upper = 0.0;
  error_code = proxy->get_tabf_alpha_boundaries(&lower, &upper);
  if (error_code != COLVARS_OK) return proxy->report_tabf_error(error_code, "get_tabf_alpha_boundaries");

  // Supply intrinsic defaults before the parent colvar parses grid parameters.
  init_scalar_boundaries(lower, upper);
  return COLVARS_OK;
}


int colvar::tabf_alpha::init_tabf(colvar &cv)
{
  // Grid parsing has already inherited the intrinsic boundaries unless overridden.
  // Even a matching user value (within tolerance) must use the exact back-end limit.
  int error_code = COLVARS_OK;
  if (cv.is_enabled(f_cv_lower_boundary) &&
      cvm::fabs(cv.lower_boundary.real_value - lower_boundary.real_value) > colvar_boundaries_tol * cv.width) {
    error_code |= cvmodule->error(
      "Error: tabfAlpha lowerBoundary " + cvm::to_str(cv.lower_boundary) +
      " conflicts with back-end TABF alpha minimum " + cvm::to_str(lower_boundary) + ".\n", COLVARS_INPUT_ERROR);
  }
  if (cv.is_enabled(f_cv_upper_boundary) &&
      cvm::fabs(cv.upper_boundary.real_value - upper_boundary.real_value) > colvar_boundaries_tol * cv.width) {
    error_code |= cvmodule->error(
      "Error: tabfAlpha upperBoundary " + cvm::to_str(cv.upper_boundary) +
      " conflicts with back-end TABF alpha maximum " + cvm::to_str(upper_boundary) + ".\n", COLVARS_INPUT_ERROR);
  }
  if (error_code != COLVARS_OK) return error_code;
  cv.lower_boundary = lower_boundary;
  cv.upper_boundary = upper_boundary;
  cv.enable(f_cv_external);

  int const factor = cv.get_time_step_factor();
  if (factor != 1) {
    return cvmodule->error("Error: timeStepFactor > 1 is not yet supported for TABF alpha variables.");
  }
  error_code = cvmodule->proxy->request_tabf_energy_freq(factor);
  if (error_code != COLVARS_OK) return cvmodule->proxy->report_tabf_error(error_code, "request_tabf_energy_freq");
  return COLVARS_OK;
}


void colvar::tabf_alpha::calc_value()
{
  // By default, follow external parameter
  // This might get overwritten by driving extended dynamics
  auto *proxy = cvmodule->proxy;
  int error_code = proxy->get_tabf_alpha(&x.real_value);
  if (error_code != COLVARS_OK) {
    proxy->report_tabf_error(error_code, "get_tabf_alpha");
    return;
  }
  error_code = proxy->get_dE_dtabf_alpha(&ft.real_value);
  if (error_code != COLVARS_OK) {
    proxy->report_tabf_error(error_code, "get_dE_dtabf_alpha");
    return;
  }
  ft.real_value *= -1.0; // Convert energy derivative to force
}


void colvar::tabf_alpha::calc_force_invgrads()
{
  // All the work is done in calc_value()
}


void colvar::tabf_alpha::calc_Jacobian_derivative()
{
  jd = 0.0;
}


void colvar::tabf_alpha::apply_force(colvarvalue const & /* force */)
{
  // Forces, if any, are applied in colvar::update_extended_Lagrangian()
}


int colvar::tabf_alpha::set_external_value(bool now)
{
  cvmodule->proxy->set_tabf_alpha(x.real_value);
  if (now) {
    return cvmodule->proxy->flush_tabf_alpha();
  }
  return COLVARS_OK;
}
