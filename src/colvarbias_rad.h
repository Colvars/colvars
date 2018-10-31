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

  /// \brief Restraint centers
  std::vector<colvarvalue> colvar_centers;

  /// \brief Restraint original centers
  std::vector<colvarvalue> colvar_orig_centers;

  /// \brief Restraint widths
  std::vector<cvm::real> colvar_widths_c;

  /// \brief Errors associated with the restraint centers
  std::vector<cvm::real> colvar_centers_errors;

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

  /// experimental files for time dependent observables (e.g. DEER)

  std::vector<std::string> time_files;

  /// experimental times for time dependent variables

  std::vector<colvarvalue> colvar_times;

  /// experimental values for time dependent variables

  std::vector<colvarvalue> colvar_expval;

  /// Specify type of colvar for parameters optimization

  std::vector<std::string> colvar_types;

  /// number of variable of a certain type

  std::vector<int> numtypes;

  cvm::matrix2d<int> whichtypes;

  /// optimize paramters

  bool opt_params;

  /// write output

  bool rad_out;

  std::ofstream radoutfile;
  std::ofstream radparfile;
  std::string rad_out_file;
  std::string rad_par_file;
  int rad_out_freq;

  /// equilibrations steps

  int equil_steps;

  /// whether set starting parameters from initial equilibration

  bool set_params_equil;

  /// mdepth of deer

  std::vector<cvm::real> mdepth_deer;

  /// backgroud parameter of deer

  std::vector<cvm::real> alpha_deer;

  /// maximum and minimum values of error scale

  cvm::real colvar_maxerror_scale;

  /// average experimental error

  cvm::real colvar_cum_error;

  /// average scale width

  cvm::real colvar_cum_uscale;

  /// unique error scale

  cvm::real colvar_errors_scale;

  /// local deviation from experiment

  std::vector<colvarvalue> colvar_deviation;

  /// total size of collective variables

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
