// src/colvarcomp_harmonicforceconstant.h

#ifndef COLVARCOMP_HARMONICFORCECONSTANT_H
#define COLVARCOMP_HARMONICFORCECONSTANT_H

#include "colvar.h"
#include "colvarcomp.h"

// Forward declaration to avoid circular dependency with the bias class.
class colvarbias_restraint;

/// \brief A fictitious scalar coordinate (lambda in [0,1]) intended for extended-Lagrangian dynamics.
///
/// This CVC does not depend on atomic coordinates and provides no gradients.
/// It is meant to be used as an extended-Lagrangian CV whose value is integrated
/// by colvar::update_extended_Lagrangian().
/// Forces from biases are applied via colvar::add_bias_force() on the parent colvar.
class cvc_harmonicforceconstant : public colvar::cvc {
public:
  cvc_harmonicforceconstant();
  virtual int init(std::string const &conf) override;
  virtual ~cvc_harmonicforceconstant() {}

  virtual void calc_value();
  virtual void calc_gradients();
  
  /// No total/system force is computed for this fictitious coordinate.
  virtual void calc_force_invgrads(); 

  /// The force applied to this CVC is handled by the parent colvar's extended Lagrangian dynamics.
  virtual void apply_force(colvarvalue const &force);
  
  /// The Jacobian derivative for this 1D fictitious coordinate is zero.
  virtual void calc_Jacobian_derivative();
  
  /// Store the value locally; do NOT talk to the MD engine.
  virtual void set_value(colvarvalue const &new_value, bool now=false) override;
  
protected:
};

#endif
