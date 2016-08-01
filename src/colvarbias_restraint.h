// -*- c++ -*-

#ifndef COLVARBIAS_RESTRAINT_H
#define COLVARBIAS_RESTRAINT_H

#include "colvarbias.h"

/// \brief Bias restraint, optionally moving towards a target
/// (implementation of \link colvarbias \endlink)
class colvarbias_restraint : public colvarbias {

public:

  /// Retrieve colvar values and calculate their biasing forces
  virtual int update();

  // TODO the following can be supplanted by a new call to init()
  /// Load new configuration - force constant and/or centers only
  virtual void change_configuration(std::string const &conf);

  /// Calculate change in energy from using alternate configuration
  virtual cvm::real energy_difference(std::string const &conf);

  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &conf);
  // virtual std::ostream & write_state_data(std::ostream &os);
  // virtual std::istream & read_state_data(std::istream &os);
  virtual std::ostream & write_state(std::ostream &os);
  virtual std::istream & read_state(std::istream &is);

  virtual std::ostream & write_traj_label(std::ostream &os);
  virtual std::ostream & write_traj(std::ostream &os);

  /// \brief Constructor
  colvarbias_restraint(char const *key);

  virtual int init(std::string const &conf);
  virtual ~colvarbias_restraint();


protected:

  /// \brief Potential function
  virtual cvm::real restraint_potential(cvm::real k, colvar const *x,
                                        colvarvalue const &xcenter) const = 0;

  /// \brief Force function
  virtual colvarvalue restraint_force(cvm::real k, colvar const *x,
                                      colvarvalue const &xcenter) const = 0;

  ///\brief Unit scaling
  virtual cvm::real restraint_convert_k(cvm::real k, cvm::real dist_measure) const = 0;

  /// \brief Restraint centers
  std::vector<colvarvalue> colvar_centers;

  /// \brief Restraint centers without wrapping or constraints applied
  std::vector<colvarvalue> colvar_centers_raw;

  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &conf);

  /// \brief Moving target?
  bool b_chg_centers;

  /// \brief New restraint centers
  std::vector<colvarvalue> target_centers;

  /// \brief Amplitude of the restraint centers' increment at each step
  /// (or stage) towards the new values (calculated from target_nsteps)
  std::vector<colvarvalue> centers_incr;

  /// Whether to write the current restraint centers to the trajectory file
  bool b_output_centers;

  /// Whether to write the current accumulated work to the trajectory file
  bool b_output_acc_work;

  /// \brief Accumulated work
  cvm::real acc_work;

  /// \brief Restraint force constant
  cvm::real force_k;

  /// \brief Changing force constant?
  bool b_chg_force_k;

  /// \brief Restraint force constant (target value)
  cvm::real target_force_k;

  /// \brief Restraint force constant (starting value)
  cvm::real starting_force_k;

  /// \brief Lambda-schedule for custom varying force constant
  std::vector<cvm::real> lambda_schedule;

  /// \brief Exponent for varying the force constant
  cvm::real force_k_exp;

  /// \brief Intermediate quantity to compute the restraint free energy
  /// (in TI, would be the accumulating FE derivative)
  cvm::real restraint_FE;


  /// \brief Equilibration steps for restraint FE calculation through TI
  cvm::real target_equil_steps;

  /// \brief Number of stages over which to perform the change
  /// If zero, perform a continuous change
  int target_nstages;

  /// \brief Number of current stage of the perturbation
  int stage;

  /// \brief Number of steps required to reach the target force constant
  /// or restraint centers
  long target_nsteps;
};


/// \brief Harmonic bias restraint
/// (implementation of \link colvarbias_restraint \endlink)
class colvarbias_restraint_harmonic : public colvarbias_restraint {

public:

  colvarbias_restraint_harmonic(char const *key);
  virtual int init(std::string const &conf);
  // no additional members, destructor not needed

protected:

  /// \brief Potential function
  virtual cvm::real restraint_potential(cvm::real k, colvar const *x,
                                        colvarvalue const &xcenter) const;

  /// \brief Force function
  virtual colvarvalue restraint_force(cvm::real k, colvar const *x,
                                      colvarvalue const &xcenter) const;

  ///\brief Unit scaling
  virtual cvm::real restraint_convert_k(cvm::real k, cvm::real dist_measure) const;

};


/// \brief Linear bias restraint
/// (implementation of \link colvarbias_restraint \endlink)
class colvarbias_restraint_linear : public colvarbias_restraint {

public:
  colvarbias_restraint_linear(char const *key);
  virtual int init(std::string const &conf);
  // no additional members, destructor not needed

protected:

  /// \brief Potential function
  virtual cvm::real restraint_potential(cvm::real k, colvar const *x,
                                        colvarvalue const &xcenter) const;

  /// \brief Force function
  virtual colvarvalue restraint_force(cvm::real k, colvar const *x,
                                      colvarvalue const &xcenter) const;

  ///\brief Unit scaling
  virtual cvm::real restraint_convert_k(cvm::real k, cvm::real dist_measure) const;

};


#endif
