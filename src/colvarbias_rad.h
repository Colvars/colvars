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
  virtual int update();
  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &state_conf);

  /// \brief Restraint centers
  std::vector<colvarvalue> colvar_centers;

  /// \brief Errors associated with the restraint centers
  std::vector<cvm::real> colvar_centers_errors;

  enum kernel_type_e {
    kt_none,
    kt_uniform,
    kt_inv_sqrt_time,
    kt_ntot
  };

  /// Specify which type of averaging is being done
  int kernel_type;

  /// Characteristic time for the kernel
  cvm::real kernel_coupling_time;

  /// Total number of samples included in the kernel averaging
  size_t kernel_num_samples;

  /// Mean deviations according to the chosen kernel times number of samples
  std::vector<colvarvalue> colvar_total_deviations;

};

#endif
