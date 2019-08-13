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

class colvarbias_reweightaMD : public colvarbias_histogram {
public:
  colvarbias_reweightaMD(char const *key);
  ~colvarbias_reweightaMD();
  virtual int init(std::string const &conf);
  virtual int update();
  virtual int write_output_files();
  void counts_to_pmf(std::vector<cvm::real>& counts) const;
  std::vector<cvm::real> compute_cumulant_expansion_factor(const std::vector<cvm::real>& dV, const std::vector<cvm::real>& dV_square, const std::vector<cvm::real>& count, cvm::real beta) const;
protected:
  /// Current accelMD factor is the from previous frame
  std::vector<int> previous_bin;
  /// Start collecting samples after N steps
  size_t start_after_steps;
  /// PMF output
  std::string out_name_pmf;

  /// Use cumulant expansion to second order?
  bool use_cumulant_expansion;
  colvar_grid_scalar* grid_count;
  colvar_grid_scalar* grid_dV;
  colvar_grid_scalar* grid_dV_square;
};

#endif
