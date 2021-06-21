// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARBIAS_HISTOGRAM_H
#define COLVARBIAS_HISTOGRAM_H

#include <vector>
#include <list>
#include <sstream>
#include <iomanip>

#include "colvarbias.h"
#include "colvargrid.h"

/// Histogram "bias" (does as the name says)
class colvarbias_histogram : public colvarbias {

public:

  colvarbias_histogram(char const *key);
  ~colvarbias_histogram();
  virtual int init(std::string const &conf);
  virtual int update();
  virtual int write_output_files();

protected:

  /// n-dim histogram
  colvar_grid_scalar *grid;
  std::vector<int> bin;
  std::string out_name, out_name_dx;

  /// If one or more of the variables are \link colvarvalue::type_vector \endlink, treat them as arrays of this length
  size_t colvar_array_size;
  /// If colvar_array_size is larger than 1, weigh each one by this number before accumulating the histogram
  std::vector<cvm::real> weights;

  virtual std::istream & read_state_data(std::istream &is);
  virtual std::ostream & write_state_data(std::ostream &os);
};

/// Reweighted histogram for accelerated molecular dynamics
class colvarbias_reweightaMD : public colvarbias_histogram {
public:
  colvarbias_reweightaMD(char const *key);
  ~colvarbias_reweightaMD();
  virtual int init(std::string const &conf);
  virtual int update();
  virtual int write_output_files();
  void hist_to_pmf(
    colvar_grid_scalar* hist,
    const colvar_grid_scalar* hist_count) const;
  void compute_cumulant_expansion_factor(
    const colvar_grid_scalar* hist_dV,
    const colvar_grid_scalar* hist_dV_square,
    const colvar_grid_scalar* hist_count,
    colvar_grid_scalar* cumulant_expansion_factor) const;
  virtual int write_exponential_reweighted_pmf(
    const std::string& output_prefix, bool append = false);
  virtual int write_cumulant_expansion_pmf(
    const std::string& output_prefix, bool append = false);
  virtual int write_count(const std::string& output_name, bool append = false);
protected:
  /// Current accelMD factor is the from previous frame
  std::vector<int> previous_bin;
  /// Start collecting samples after N steps
  size_t start_after_steps;

  /// Use cumulant expansion to second order?
  bool b_use_cumulant_expansion;
  colvar_grid_scalar* grid_count;
  colvar_grid_scalar* grid_dV;
  colvar_grid_scalar* grid_dV_square;

  /// Number of timesteps between recording data in history files (if non-zero)
  size_t history_freq;
  bool b_history_files;

  /// Write gradients of the PMF?
  bool b_write_gradients;

  /// save and restore
  virtual std::istream & read_state_data(std::istream &is) override;
  virtual std::ostream & write_state_data(std::ostream &os) override;
private:
  /// temporary grids for evaluating PMFs
  colvar_grid_scalar  *pmf_grid_exp_avg;
  colvar_grid_scalar  *pmf_grid_cumulant;
  colvar_grid_gradient *grad_grid_exp_avg;
  colvar_grid_gradient *grad_grid_cumulant;
};

#endif
