// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARBIAS_ABMD_H
#define COLVARBIAS_ABMD_H

#include "colvarbias_restraint.h"


/// \brief Adiabatic Bias MD
class colvarbias_abmd
  : public virtual colvarbias,
    public virtual colvarbias_ti
{
public:

  colvarbias_abmd(colvarmodule *cvmodule_in, char const *key);
  ~colvarbias_abmd() = default;

  int init(std::string const &conf) override;
  int update() override;
  std::string const get_state_params() const  override;
  int set_state_params(std::string const &conf) override;
  std::ostream & write_traj_label(std::ostream &os) override;
  std::ostream & write_traj(std::ostream &os) override;

protected:

  /// \brief Location of the moving wall
  cvm::real ref_val = 0.;

  /// \brief Has ref_val already been set?
  bool ref_initialized = false;

  /// \brief Value of the reference where it stops moving
  cvm::real stopping_val = 0.;

  /// \brief Is the target moving down?
  bool decreasing = false;

  /// \brief Restraint force constant
  cvm::real k = 0.;
};


#endif
