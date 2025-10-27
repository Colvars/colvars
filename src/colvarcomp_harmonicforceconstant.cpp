// src/colvarcomp_harmonicforceconstant.cpp

#include "colvarcomp_harmonicforceconstant.h"
#include "colvarbias_restraint.h"
#include "colvarmodule.h"

cvc_harmonicforceconstant::cvc_harmonicforceconstant()
    : cvc() // Call base class constructor
{
  set_function_type("harmonicForceConstant"); // Set the type name for this component
  // This is a fictitious coordinate, so it does not have gradients with respect to atomic positions.
  provide(f_cvc_explicit_gradient, false);
  provide(f_cvc_gradient, false);
  provide(f_cvc_collect_atom_ids, false);
  // It does, however, provide a "total force" which is the thermodynamic force F_lambda = -dU/d_lambda.
  provide(f_cvc_inv_gradient);
  // The Jacobian derivative is zero for a 1D fictitious coordinate.
  provide(f_cvc_Jacobian);
  
  harmonic_bias = NULL;
  is_linked = false;
  k_exponent = 1.0; // Default to linear scaling

  // The colvar object controls the extended Lagrangian dynamics via the f_cv_external flag.

  // This CVC is a scalar defined on the interval [0, 1].
  init_scalar_boundaries(0.0, 1.0); 
  
  x.type(colvarvalue::type_scalar);
  x.real_value = 0.0; // Initial value, corresponding to k=0.

  cvm::log("Initializing a harmonicForceConstant component.\n");
}

int cvc_harmonicforceconstant::init(std::string const &conf)
{
  cvc::init(conf);
  if (!get_keyval(conf, "harmonicName", harmonic_bias_name, std::string(""))) {
    cvm::error("Error: Missing required parameter harmonicName for harmonicForceConstant component.");
    return COLVARS_INPUT_ERROR;
  }
  
  // Parse the optional kExponent parameter
  if (get_keyval(conf, "kExponent", k_exponent, k_exponent)) {
    if (k_exponent <= 0.0) {
        cvm::error("Error: kExponent must be positive for harmonicForceConstant component.", COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
    }
    cvm::log("Using exponent kExponent = " + cvm::to_str(k_exponent) + " for force constant scaling.\n");
  }
  
  return COLVARS_OK;
}

int cvc_harmonicforceconstant::link_bias(colvarmodule *cvm, colvar *cv)
{
  if (is_linked) return COLVARS_OK; // Already linked, do nothing.

  colvarbias *bias = cvm->bias_by_name(harmonic_bias_name);
  if (!bias) {
    cvm->error("Error: Cannot find harmonic bias named '" + harmonic_bias_name + "' for harmonicForceConstant component.");
    return COLVARS_INPUT_ERROR;
  }
  
  // Attempt to cast the generic bias to a restraint.
  harmonic_bias = dynamic_cast<colvarbias_restraint *>(bias);
  if (!harmonic_bias) {
    cvm->error("Error: Bias '" + harmonic_bias_name + "' is not a harmonic restraint.");
    return COLVARS_INPUT_ERROR;
  }
  
  // Register this CV with the harmonic bias to enable dynamic control.
  // The 'cv' parameter is the parent colvar that owns this CVC.
  harmonic_bias->set_dynamic_k_cv(cv);
  is_linked = true;
  cvm::log("Successfully linked harmonicForceConstant component to harmonic bias '" + harmonic_bias_name + "'.\n");

  return COLVARS_OK;
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

void cvc_harmonicforceconstant::calc_force_invgrads()
{
  // This function computes the "total force" on the fictitious coordinate, which is F_lambda = -dU/d_lambda.
  // This force is calculated and cached within the controlled harmonic bias.
  if (is_linked && harmonic_bias) {
    // Retrieve F_lambda calculated in the previous timestep by the bias
    // and store it in this CVC's 'ft' (total_force) member.
    ft.real_value = harmonic_bias->get_k_derivative();
  } else {
    ft.real_value = 0.0;
  }
}

void cvc_harmonicforceconstant::calc_Jacobian_derivative() {
  // The Jacobian derivative (metric correction) for a 1D fictitious coordinate is zero.
  jd.type(colvarvalue::type_scalar);
  jd.real_value = 0.0;
}

void cvc_harmonicforceconstant::apply_force(colvarvalue const &force)
{
  // This CVC does not interact with any atoms directly.
  // Forces from biases (ABF, Metadynamics) are applied to the extended coordinate
  // in colvar::update_forces_energy(). The thermodynamic force from the restraint
  // is reported via calc_force_invgrads().
}