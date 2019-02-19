// -*- c++ -*-

#ifndef COLVARBIAS_RAD_H
#define COLVARBIAS_RAD_H

#include "colvarbias.h"


/// \brief Restrained Average Dynamics (RAD)
/// Marinelli and Fiorin, Structure, 2018

class colvarbias_rad : public colvarbias {

public:

  colvarbias_rad(char const *key);
  virtual int init(std::string const &conf);
  virtual int init_centers(std::string const &conf);
  virtual int clear();
  virtual int update();
  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &state_conf);
  virtual int setup_output();
  virtual int write_traj_files();

  /// \brief Restraint centers
  std::vector<colvarvalue> colvar_centers;

  /// \brief Restraint original centers
  std::vector<colvarvalue> colvar_exp_centers;

  /// \brief Scale factors for the restraint
  std::vector<cvm::real> colvar_widths_c;

  /// Average of colvar_widths_c
  cvm::real colvar_cum_uscale;

  /// \brief Errors associated with the restraint centers
  std::vector<cvm::real> colvar_centers_errors;

  /// \brief value of the experimental parameters associated to each CV

  std::vector<cvm::vector1d<cvm::real> > val_params;

  /// Types of averaging
  enum kernel_type_e {
    kt_none,
    kt_uniform,
    kt_inv_sqrt_time,
    kt_ntot
  };

  /// Specify which type of averaging is being done
  kernel_type_e kernel_type;

  /// Characteristic time for the kernel
  cvm::real kernel_coupling_time;

  /// Characteristic time for parameters optimization
  cvm::real params_coupling_time;

  /// frequency at which to write output
  int rad_out_freq;

  /// write progressive output of RAD convergence
  std::ostream *rad_out_os;

  /// Name of output file
  inline std::string rad_out_file_name() const
  {
    return cvm::output_prefix() + "." + this->name + ".traj";
  }

  /// Types of parameters optimization of the associated variables

  enum opt_type_e {
    opt_none,
    opt_lambda,
    opt_chiquare
  };

  /// Number of optimization steps for parameters optimization using chisquare minimization
  size_t colvar_chisquare_opt_steps;

  size_t colvar_rad_steps; // XXX just for checking to be deleted after debugging

  /// Mean total deviation respect to the experimental values:
  /// useful to optimize the parameters by gradient based
  /// minimization of the chisquare function
  std::vector<colvarvalue> colvar_total_chideviations;

  /// lambda based parameters optimization
  bool lambda_par_opt;

  /// parameters optimization based on chi square minimization
  bool chi_square_par_opt;

  /// write progressive output of CVC parameter optimization by RAD
  std::ostream *rad_param_os;

  /// Name of parameter optimization output file
  inline std::string rad_param_file_name() const
  {
    return cvm::output_prefix() + "." + this->name + ".params.traj";
  }

  /// equilibration steps
  int equil_steps;

  /// average experimental error
  cvm::real colvar_cum_error;

  /// unique error scale
  cvm::real colvar_errors_scale;

  /// local deviation from experiment
  std::vector<colvarvalue> colvar_deviation;

  /// average deviation from experiment
  cvm::real colvar_aver_deviation;

  /// total dimensionality of all involved collective variables
  int colvar_size_tot;

  /// update unique error scale to fix chi to one
  bool fix_chi_square_one;

  /// update individual error scales to fix chi to one
  bool fix_chi_square_multi;

  /// use norm 1 to calculate deviation from experiment
  bool use_norm_1;

  /// Mean of the deviation according to the chosen kernel times number of samples
  std::vector<colvarvalue> colvar_total_deviations;

};

#endif
