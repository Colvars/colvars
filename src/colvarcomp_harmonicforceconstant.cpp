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

void cvc_harmonicforceconstant::calc_force_invgrads()
{
  // Default to zero force
  ft.real_value = 0.0;
  // Find the harmonic bias this CVC is linked to
  colvarbias *bias = cvm::main()->bias_by_name(this->harmonic_bias_name);
  if (!bias) {
    // This might happen during initialization, it's not an error yet.
    return;
  }
  colvarbias_restraint_harmonic* h_bias = dynamic_cast<colvarbias_restraint_harmonic*>(bias);
  if (!h_bias) {
    cvm::error("Error: Bias '" + this->harmonic_bias_name + "' is not a harmonic restraint.", COLVARS_INPUT_ERROR);
    return;
  }
  // Get the maximum force constant from the bias
  // NOTE: This requires access to a protected member. A cleaner solution would be
  // to add a public getter `get_force_k()` to `colvarbias_restraint_k`.
  // For now, we assume we can get it. Let's add that getter.
  // In `colvarbias_restraint.h`, in `colvarbias_restraint_k`, add:
  //   cvm::real get_force_k() const { return force_k; }
  
  cvm::real k_max = h_bias->get_force_k(); // Assuming you add this getter
  cvm::real raw_lambda = parent->value().real_value; // Use the parent colvar's value
  // Calculate thermodynamic force F_lambda = -dU/d_lambda
  cvm::real dU_d_k_eff = h_bias->get_dU_d_k_eff();
  cvm::real d_k_eff_d_lambda;
  if (raw_lambda == 0.0) {
      if (k_exponent < 1.0) d_k_eff_d_lambda = std::numeric_limits<cvm::real>::infinity();
      else if (k_exponent == 1.0) d_k_eff_d_lambda = k_max;
      else d_k_eff_d_lambda = 0.0;
  } else {
      d_k_eff_d_lambda = k_max * k_exponent * cvm::pow(raw_lambda, k_exponent - 1.0);
  }
  if (!std::isfinite(d_k_eff_d_lambda)) {
      cvm::log("Warning: Derivative factor for k is non-finite in CVC. Setting thermodynamic force to 0.");
      d_k_eff_d_lambda = 0.0;
  }
  // F_lambda = -dU/d_lambda = - (dU/dk_eff) * (dk_eff/d_lambda)
  ft.real_value = -1.0 * dU_d_k_eff * d_k_eff_d_lambda;
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