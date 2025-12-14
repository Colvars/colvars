// src/colvarcomp_harmonicforceconstant.cpp

#include "colvarcomp_harmonicforceconstant.h"
#include "colvarmodule.h"

cvc_harmonicforceconstant::cvc_harmonicforceconstant()
    : cvc() // Call base class constructor
{
  set_function_type("harmonicForceConstant"); // Set the type name for this component
  // This is a fictitious coordinate, so it does not have gradients with respect to atomic positions.
  provide(f_cvc_explicit_gradient, false);
  provide(f_cvc_gradient, false);
  provide(f_cvc_collect_atom_ids, false);
  // This fictitious coordinate does not support total-force calculation
  // through inverse gradients / Jacobian machinery.
  provide(f_cvc_inv_gradient, false);
  provide(f_cvc_Jacobian, false);

  // The colvar object controls the extended Lagrangian dynamics via the f_cv_external flag.

  // This CVC is a scalar defined on the interval [0, 1].
  init_scalar_boundaries(0.0, 1.0); 
  
  x.type(colvarvalue::type_scalar);
  x.real_value = 0.0; // Initial value, corresponding to k=0.

  cvm::log("Initializing a harmonicForceConstant component.\n");
}

int cvc_harmonicforceconstant::init(std::string const &conf)
{
  // No parameters: behave like alchLambda-style external coordinate
  return cvc::init(conf);
}

void cvc_harmonicforceconstant::calc_force_invgrads()
{
  // Not implemented: this coordinate has no inverse gradients.
}

void cvc_harmonicforceconstant::calc_value()
{
  // The value of this CV is managed by the parent colvar's extended Lagrangian dynamics.
  // This function does nothing.
}

void cvc_harmonicforceconstant::calc_gradients()
{
  // This is a fictitious coordinate and has no gradients with respect to atomic positions.
  // This function does nothing.
}

void cvc_harmonicforceconstant::calc_Jacobian_derivative() {
  // Not implemented
  cvm::error("Error: Jacobian derivative is not implemented for harmonicForceConstant.\n",
             COLVARS_NOT_IMPLEMENTED);
}

void cvc_harmonicforceconstant::apply_force(colvarvalue const &force)
{
  // This CVC does not interact with any atoms directly.
  // Forces from biases (ABF, Metadynamics) are applied to the extended coordinate
  // in colvar::update_forces_energy(). The thermodynamic force from the restraint
  // is reported via calc_force_invgrads().
}

void cvc_harmonicforceconstant::set_value(colvarvalue const &new_value, bool now)
{
  x = new_value;
}