// src/colvarcomp_harmonicforceconstant.h

#ifndef COLVARCOMP_HARMONICFORCECONSTANT_H
#define COLVARCOMP_HARMONICFORCECONSTANT_H

#include "colvar.h"
#include "colvarcomp.h"

// Forward declaration to avoid circular dependency with the bias class.
class colvarbias_restraint;

/// \brief A fictitious coordinate representing the normalized force constant of a harmonic restraint.
///
/// This CVC allows a harmonic restraint's force constant (k) to be treated as a
/// dynamic variable, lambda, evolving between 0 and 1. The potential energy of the
/// restraint is defined as U(x, lambda) = 0.5 * (k_max * lambda^n) * (x-x_0)^2.
/// This CVC reports the thermodynamic force F_lambda = -dU/d_lambda, enabling
/// enhanced sampling methods (like ABF or Metadynamics) to reconstruct the free energy
/// profile associated with growing or disappearing the restraint.
class cvc_harmonicforceconstant : public colvar::cvc {
public:
  cvc_harmonicforceconstant();
  virtual int init(std::string const &conf);
  virtual ~cvc_harmonicforceconstant() {}

  /// Called after all colvars and biases are initialized to link to the target restraint.
  virtual int link_bias(colvarmodule *cvm, colvar *cv);

  virtual void calc_value();
  virtual void calc_gradients();
  
  /// This CVC reports the thermodynamic force F_lambda, so this function is implemented.
  virtual void calc_force_invgrads(); 

  /// The force applied to this CVC is handled by the parent colvar's extended Lagrangian dynamics.
  virtual void apply_force(colvarvalue const &force);
  
  /// The Jacobian derivative for this 1D fictitious coordinate is zero.
  virtual void calc_Jacobian_derivative();
  
  /// Returns the exponent used for scaling the force constant.
  cvm::real get_k_exponent() const { return k_exponent; }

protected:
  /// The name of the harmonic bias to be controlled.
  std::string harmonic_bias_name;
  /// A pointer to the instance of the harmonic bias.
  colvarbias_restraint *harmonic_bias;
  
  /// Exponent 'n' for the force constant scaling, k = k_max * lambda^n.
  cvm::real k_exponent = 1.0;
  
  /// Flag to ensure the link_bias function is executed only once.
  bool is_linked;
};

#endif