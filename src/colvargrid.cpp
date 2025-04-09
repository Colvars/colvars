// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <ctime>

#include "colvar.h"
#include "colvargrid.h"
#include "colvargrid_def.h"
#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvarvalue.h"

#include <fstream>
#include <sstream>

#include <algorithm>
#include <iostream>

// Helper function to print vector<int>
std::string vec_to_string(const std::vector<int> &vec)
{
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    oss << vec[i];
    if (i < vec.size() - 1)
      oss << ", ";
  }
  oss << "]";
  return oss.str();
};

colvar_grid_count::colvar_grid_count() : colvar_grid<size_t>() { mult = 1; }

colvar_grid_count::colvar_grid_count(std::vector<colvar *> &colvars, std::string config)
    : colvar_grid<size_t>(colvars, 0, 1, false, nullptr, config)
{
}

colvar_grid_count::colvar_grid_count(std::vector<colvar *> &colvars,
                                     std::shared_ptr<const colvar_grid_params> params)
    : colvar_grid<size_t>(colvars, 0, 1, false, params)
{
}

colvar_grid_count::colvar_grid_count(std::string &filename)
{
  // multiplicity is 1 by definition
  init_from_file(filename, 1);
}

std::string colvar_grid_count::get_state_params() const
{
  return colvar_grid<size_t>::get_state_params();
}

int colvar_grid_count::parse_params(std::string const &conf,
                                    colvarparse::Parse_Mode const parse_mode)
{
  return colvar_grid<size_t>::parse_params(conf, parse_mode);
}

std::istream &colvar_grid_count::read_restart(std::istream &is)
{
  return colvar_grid<size_t>::read_restart(is);
}

cvm::memory_stream &colvar_grid_count::read_restart(cvm::memory_stream &is)
{
  return colvar_grid<size_t>::read_restart(is);
}

std::ostream &colvar_grid_count::write_restart(std::ostream &os)
{
  return colvar_grid<size_t>::write_restart(os);
}

cvm::memory_stream &colvar_grid_count::write_restart(cvm::memory_stream &os)
{
  return colvar_grid<size_t>::write_restart(os);
}

std::istream &colvar_grid_count::read_raw(std::istream &is)
{
  return colvar_grid<size_t>::read_raw(is);
}

cvm::memory_stream &colvar_grid_count::read_raw(cvm::memory_stream &is)
{
  return colvar_grid<size_t>::read_raw(is);
}

std::ostream &colvar_grid_count::write_raw(std::ostream &os, size_t const buf_size) const
{
  return colvar_grid<size_t>::write_raw(os, buf_size);
}

cvm::memory_stream &colvar_grid_count::write_raw(cvm::memory_stream &os,
                                                 size_t const buf_size) const
{
  return colvar_grid<size_t>::write_raw(os, buf_size);
}

std::istream &colvar_grid_count::read_multicol(std::istream &is, bool add)
{
  return colvar_grid<size_t>::read_multicol(is, add);
}

int colvar_grid_count::read_multicol(std::string const &filename, std::string description, bool add)
{
  return colvar_grid<size_t>::read_multicol(filename, description, add);
}

std::ostream &colvar_grid_count::write_multicol(std::ostream &os) const
{
  return colvar_grid<size_t>::write_multicol(os);
}

int colvar_grid_count::write_multicol(std::string const &filename, std::string description) const
{
  return colvar_grid<size_t>::write_multicol(filename, description);
}

std::ostream &colvar_grid_count::write_opendx(std::ostream &os) const
{
  return colvar_grid<size_t>::write_opendx(os);
}

int colvar_grid_count::write_opendx(std::string const &filename, std::string description) const
{
  return colvar_grid<size_t>::write_opendx(filename, description);
}


colvar_grid_scalar::colvar_grid_scalar() : colvar_grid<cvm::real>(), samples(NULL) {}

colvar_grid_scalar::colvar_grid_scalar(colvar_grid_scalar const &g)
    : colvar_grid<cvm::real>(g), samples(NULL)
{
}

colvar_grid_scalar::colvar_grid_scalar(std::vector<colvar *> &colvars,
                                       std::shared_ptr<const colvar_grid_params> params,
                                       bool add_extra_bin)
    : colvar_grid<cvm::real>(colvars, 0.0, 1, add_extra_bin, params), samples(NULL)
{
}

colvar_grid_scalar::~colvar_grid_scalar() {}

std::string colvar_grid_scalar::get_state_params() const
{
  return colvar_grid<cvm::real>::get_state_params();
}

int colvar_grid_scalar::parse_params(std::string const &conf,
                                     colvarparse::Parse_Mode const parse_mode)
{
  return colvar_grid<cvm::real>::parse_params(conf, parse_mode);
}

std::istream &colvar_grid_scalar::read_restart(std::istream &is)
{
  return colvar_grid<cvm::real>::read_restart(is);
}

cvm::memory_stream &colvar_grid_scalar::read_restart(cvm::memory_stream &is)
{
  return colvar_grid<cvm::real>::read_restart(is);
}

std::ostream &colvar_grid_scalar::write_restart(std::ostream &os)
{
  return colvar_grid<cvm::real>::write_restart(os);
}

cvm::memory_stream &colvar_grid_scalar::write_restart(cvm::memory_stream &os)
{
  return colvar_grid<cvm::real>::write_restart(os);
}

std::istream &colvar_grid_scalar::read_raw(std::istream &is)
{
  return colvar_grid<cvm::real>::read_raw(is);
}

cvm::memory_stream &colvar_grid_scalar::read_raw(cvm::memory_stream &is)
{
  return colvar_grid<cvm::real>::read_raw(is);
}

std::ostream &colvar_grid_scalar::write_raw(std::ostream &os, size_t const buf_size) const
{
  return colvar_grid<cvm::real>::write_raw(os, buf_size);
}

cvm::memory_stream &colvar_grid_scalar::write_raw(cvm::memory_stream &os,
                                                  size_t const buf_size) const
{
  return colvar_grid<cvm::real>::write_raw(os, buf_size);
}

std::istream &colvar_grid_scalar::read_multicol(std::istream &is, bool add)
{
  return colvar_grid<cvm::real>::read_multicol(is, add);
}

int colvar_grid_scalar::read_multicol(std::string const &filename, std::string description,
                                      bool add)
{
  return colvar_grid<cvm::real>::read_multicol(filename, description, add);
}

std::ostream &colvar_grid_scalar::write_multicol(std::ostream &os) const
{
  return colvar_grid<cvm::real>::write_multicol(os);
}

int colvar_grid_scalar::write_multicol(std::string const &filename, std::string description) const
{
  return colvar_grid<cvm::real>::write_multicol(filename, description);
}

std::ostream &colvar_grid_scalar::write_opendx(std::ostream &os) const
{
  return colvar_grid<cvm::real>::write_opendx(os);
}

int colvar_grid_scalar::write_opendx(std::string const &filename, std::string description) const
{
  return colvar_grid<cvm::real>::write_opendx(filename, description);
}


cvm::real colvar_grid_scalar::maximum_value() const
{
  cvm::real max = data[0];
  for (size_t i = 0; i < nt; i++) {
    if (data[i] > max)
      max = data[i];
  }
  return max;
}


cvm::real colvar_grid_scalar::minimum_value() const
{
  cvm::real min = data[0];
  for (size_t i = 0; i < nt; i++) {
    if (data[i] < min)
      min = data[i];
  }
  return min;
}

cvm::real colvar_grid_scalar::minimum_pos_value() const
{
  cvm::real minpos = data[0];
  size_t i;
  for (i = 0; i < nt; i++) {
    if (data[i] > 0) {
      minpos = data[i];
      break;
    }
  }
  for (i = 0; i < nt; i++) {
    if (data[i] > 0 && data[i] < minpos)
      minpos = data[i];
  }
  return minpos;
}

cvm::real colvar_grid_scalar::integral() const
{
  cvm::real sum = 0.0;
  for (size_t i = 0; i < nt; i++) {
    sum += data[i];
  }
  cvm::real bin_volume = 1.0;
  for (size_t id = 0; id < widths.size(); id++) {
    bin_volume *= widths[id];
  }
  return bin_volume * sum;
}


cvm::real colvar_grid_scalar::entropy() const
{
  cvm::real sum = 0.0;
  for (size_t i = 0; i < nt; i++) {
    if (data[i] > 0) {
      sum += -1.0 * data[i] * cvm::logn(data[i]);
    }
  }
  cvm::real bin_volume = 1.0;
  for (size_t id = 0; id < widths.size(); id++) {
    bin_volume *= widths[id];
  }
  return bin_volume * sum;
}

/// \brief Return the RMSD between this grid's data and another one
/// Grids must have the same dimensions.
cvm::real colvar_grid_scalar::grid_rmsd(colvar_grid_scalar const &other_grid) const
{
  if (other_grid.data.size() != this->data.size()) {
    cvm::error("Error: trying to subtract two grids with "
               "different size.\n");
    return -1.;
  }

  cvm::real sum2 = 0.0;

  if (samples && other_grid.samples) {
    for (size_t i = 0; i < data.size(); i++) {
      size_t n = samples->get_value(i);
      cvm::real us = n ? data[i] / n : 0.0;
      n = other_grid.samples->get_value(i);
      cvm::real them = n ? other_grid.data[i] / n : 0.0;
      cvm::real d = us - them;
      sum2 += d * d;
    }
  } else {
    for (size_t i = 0; i < data.size(); i++) {
      cvm::real d = other_grid.data[i] - data[i];
      sum2 += d * d;
    }
  }

  return sqrt(sum2 / this->data.size());
}


colvar_grid_gradient::colvar_grid_gradient() : colvar_grid<cvm::real>(), samples(NULL) {}


// colvar_grid_gradient::colvar_grid_gradient(std::vector<colvar *> &colvars, std::string config)
//   : colvar_grid<cvm::real>(colvars, 0.0, colvars.size(), false, nullptr, config), samples(NULL)
// {}

// colvar_grid_gradient::colvar_grid_gradient(std::vector<colvar *> &colvars,
//                                            std::shared_ptr<colvar_grid_count> samples_in)
//   : colvar_grid<cvm::real>(colvars, 0.0, colvars.size(), false, samples_in), samples(samples_in)
// {
//   if (samples_in)
//     samples_in->has_parent_data = true;
// }

colvar_grid_gradient::colvar_grid_gradient(std::vector<colvar *> &colvars,
                                           std::shared_ptr<colvar_grid_count> samples_in,
                                           std::shared_ptr<const colvar_grid_params> params,
                                           std::string config)
    : colvar_grid<cvm::real>(colvars, 0.0, colvars.size(), false, params, config),
      samples(samples_in)
{
  if (samples_in)
    samples_in->has_parent_data = true;
}

colvar_grid_gradient::colvar_grid_gradient(std::string &filename,
                                           std::shared_ptr<colvar_grid_count> samples_in)
{
  samples = samples_in;
  // Elements will be read and multiplied by sample counts through
  // colvar_grid_gradient::value_input()
  init_from_file(filename, 0); // convention: set mult to 0 for gradient
}

std::string colvar_grid_gradient::get_state_params() const
{
  return colvar_grid<cvm::real>::get_state_params();
}

int colvar_grid_gradient::parse_params(std::string const &conf,
                                       colvarparse::Parse_Mode const parse_mode)
{
  return colvar_grid<cvm::real>::parse_params(conf, parse_mode);
}

std::istream &colvar_grid_gradient::read_restart(std::istream &is)
{
  return colvar_grid<cvm::real>::read_restart(is);
}

cvm::memory_stream &colvar_grid_gradient::read_restart(cvm::memory_stream &is)
{
  return colvar_grid<cvm::real>::read_restart(is);
}

std::ostream &colvar_grid_gradient::write_restart(std::ostream &os)
{
  return colvar_grid<cvm::real>::write_restart(os);
}

cvm::memory_stream &colvar_grid_gradient::write_restart(cvm::memory_stream &os)
{
  return colvar_grid<cvm::real>::write_restart(os);
}

std::istream &colvar_grid_gradient::read_raw(std::istream &is)
{
  return colvar_grid<cvm::real>::read_raw(is);
}

cvm::memory_stream &colvar_grid_gradient::read_raw(cvm::memory_stream &is)
{
  return colvar_grid<cvm::real>::read_raw(is);
}

std::ostream &colvar_grid_gradient::write_raw(std::ostream &os, size_t const buf_size) const
{
  return colvar_grid<cvm::real>::write_raw(os, buf_size);
}

cvm::memory_stream &colvar_grid_gradient::write_raw(cvm::memory_stream &os,
                                                    size_t const buf_size) const
{
  return colvar_grid<cvm::real>::write_raw(os, buf_size);
}

std::istream &colvar_grid_gradient::read_multicol(std::istream &is, bool add)
{
  return colvar_grid<cvm::real>::read_multicol(is, add);
}

int colvar_grid_gradient::read_multicol(std::string const &filename, std::string description,
                                        bool add)
{
  return colvar_grid<cvm::real>::read_multicol(filename, description, add);
}

std::ostream &colvar_grid_gradient::write_multicol(std::ostream &os) const
{
  return colvar_grid<cvm::real>::write_multicol(os);
}

int colvar_grid_gradient::write_multicol(std::string const &filename, std::string description) const
{
  return colvar_grid<cvm::real>::write_multicol(filename, description);
}

std::ostream &colvar_grid_gradient::write_opendx(std::ostream &os) const
{
  return colvar_grid<cvm::real>::write_opendx(os);
}

int colvar_grid_gradient::write_opendx(std::string const &filename, std::string description) const
{
  return colvar_grid<cvm::real>::write_opendx(filename, description);
}


void colvar_grid_gradient::write_1D_integral(std::ostream &os)
{
  cvm::real bin, min, integral;
  std::vector<cvm::real> int_vals;

  os << "#       xi            A(xi)\n";

  if (cv.size() != 1) {
    cvm::error("Cannot write integral for multi-dimensional gradient grids.");
    return;
  }

  integral = 0.0;
  int_vals.push_back(0.0);
  min = 0.0;

  // correction for periodic colvars, so that the PMF is periodic
  cvm::real corr;
  if (periodic[0]) {
    corr = average();
  } else {
    corr = 0.0;
  }

  for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix)) {

    if (samples) {
      size_t const samples_here = samples->value(ix);
      if (samples_here)
        integral += (value(ix) / cvm::real(samples_here) - corr) * cv[0]->width;
    } else {
      integral += (value(ix) - corr) * cv[0]->width;
    }

    if (integral < min)
      min = integral;
    int_vals.push_back(integral);
  }

  bin = 0.0;
  for (int i = 0; i < nx[0]; i++, bin += 1.0) {
    os << std::setw(10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
       << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << int_vals[i] - min << "\n";
  }

  os << std::setw(10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
     << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << int_vals[nx[0]] - min
     << "\n";

  return;
}


/// \brief Return the RMSD between this grid's data and another one
/// Grids must have the same dimensions.
cvm::real colvar_grid_gradient::grid_rmsd(colvar_grid_gradient const &other_grid) const
{
  if (other_grid.multiplicity() != this->multiplicity()) {
    cvm::error("Error: trying to subtract two grids with "
               "different multiplicity.\n");
    return -1.;
  }

  if (other_grid.data.size() != this->data.size()) {
    cvm::error("Error: trying to subtract two grids with "
               "different size.\n");
    return -1.;
  }

  cvm::real sum2 = 0.0;
  std::vector<int> ix;
  size_t imult;
  for (ix = new_index(); index_ok(ix); incr(ix)) {
    for (imult = 0; imult < this->multiplicity(); imult++) {
      cvm::real d = this->value_output(ix, imult) - other_grid.value_output(ix, imult);
      sum2 += d * d;
    }
  }
  return sqrt(sum2 / this->data.size());
}


integrate_potential::integrate_potential(std::vector<colvar *> &colvars,
                                         std::shared_ptr<colvar_grid_gradient> gradients)
    : colvar_grid_scalar(colvars, gradients, true), b_smoothed(false), gradients(gradients)
{
  // parent class colvar_grid_scalar is constructed with add_extra_bin option set to true
  // hence PMF grid is wider than gradient grid if non-PBC

  if (nd > 1) {
    //TODO: restore this
    cvm::main()->cite_feature("Poisson integration of 2D/3D free energy surfaces");
    divergence.resize(computation_nt);
    div_border_supplement.resize(computation_nt);
    weights.resize(nt);
    fdiff_gradient.resize(nt * nd);

    // Compute inverse of Laplacian diagonal for Jacobi preconditioning
    // For now all code related to preconditioning is commented out
    // until a method better than Jacobi is implemented
    //     cvm::log("Preparing inverse diagonal for preconditioning...\n");
    //     inv_lap_diag.resize(nt);
    //     std::vector<cvm::real> id(nt), lap_col(nt);
    //     for (int i = 0; i < nt; i++) {
    //       if (i % (nt / 100) == 0)
    //         cvm::log(cvm::to_str(i));
    //       id[i] = 1.;
    //       atimes(id, lap_col);
    //       id[i] = 0.;
    //       inv_lap_diag[i] = 1. / lap_col[i];
    //     }
    //     cvm::log("Done.\n");
  }
}


integrate_potential::integrate_potential(std::shared_ptr<colvar_grid_gradient> gradients)
    : b_smoothed(false), gradients(gradients)
{
  nd = gradients->num_variables();
  nx = gradients->number_of_points_vec();
  widths = gradients->widths;
  periodic = gradients->periodic;

  init_computation_nx_nt();
  divergence.resize(computation_nt);
  div_border_supplement.resize(computation_nt);
  prepare_divergence_calculation();

  // Expand grid by 1 bin in non-periodic dimensions
  for (size_t i = 0; i < nd; i++) {
    if (!periodic[i])
      nx[i]++;
    // Shift the grid by half the bin width (values at edges instead of center of bins)
    lower_boundaries.push_back(gradients->lower_boundaries[i].real_value - 0.5 * widths[i]);
  }
  //TODO: ask Jérôme if this is correct --> it wasn't and now it is
  setup(nx);
  computation_grid->periodic = periodic;
  computation_grid->setup(computation_nx);
  
  if (nd > 1) {
    divergence.resize(computation_nt);
    weights.resize(nt);
    fdiff_gradient.resize(nt * nd);
  }

}


int integrate_potential::integrate(const int itmax, const cvm::real &tol, cvm::real &err,
                                   bool verbose, bool weighted)
{
  int iter = 0;
  
  if (nd == 1 && !weighted) {

    cvm::real sum = 0.0;
    cvm::real corr;
    if (periodic[0]) {
      corr = gradients->average(); // Enforce PBC by subtracting average gradient
    } else {
      corr = 0.0;
    }
    //TODO: ask Jérôme what does this do? --> integrate in dimension one = dz * value
    std::vector<int> ix;
    // Iterate over valid indices in gradient grid
    for (ix = new_index(); gradients->index_ok(ix); incr(ix)) {
      set_value(ix, sum);
      cvm::real val = gradients->value_output_smoothed(ix, b_smoothed);
      sum += (val - corr) * widths[0];
    }
    if (index_ok(ix)) {
      // This will happen if non-periodic: then PMF grid has one extra bin wrt gradient grid
      // If not, sum should be zero
      set_value(ix, sum);
    }

  } else if (nd <= 3) {
    if (weighted){
      set_weighted_div();
      laplacian_weighted<true>(divergence, data);
      for (int i = 0; i < computation_nt; i++){
        divergence[i] += div_border_supplement[i];
      }
    }

    
    nr_linbcg_sym(true, divergence, computation_grid->data, tol, itmax, iter, err);
    if (verbose)
      cvm::log("Integrated in " + cvm::to_str(iter) + " steps, error: " + cvm::to_str(err));

    // DEBUG ###########################
    auto backup = data;
    data = divergence;
    std::ofstream os("div.dat");
    write_multicol(os);
    os.close();
    data = weights;
    os.open("weights.dat");
    write_multicol(os);
    os.close();
    data = backup;
    // DEBUG 2 ###########################
    // Compute terms of the Laplacian matrix
    std::vector<cvm::real> lap_mat(nt, 0.);

    std::vector<size_t> cols = {0, 1, 2, 3, 4, 5, nt - 6, nt - 5, nt - 4, nt - 3, nt - 2, nt - 1};

    for (size_t i = 0; i < cols.size(); i++) {
      this->reset();
      data[cols[i]] = 1.;
      laplacian_weighted<true>(data, lap_mat);
      printf("Col  %3li  | ", cols[i]);
      for (size_t j = 0; j < cols.size(); j++) {
        printf(" %6.1f", lap_mat[cols[j]]);
      }
      printf("\n");
    }
    // DEBUG 2 ###########################


  } else {
    cvm::error("Cannot integrate PMF in dimension > 3\n");
  }

  return iter;
}


void integrate_potential::set_div()
{
  if (nd == 1)
    return;

  for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix)) {
      update_weighted_div_local(ix);
  }
}

void integrate_potential::set_weighted_div()
  {
    sum_count = 0;
    std::vector<int> max_position;
    std::vector<int> min_position;
    sorted_counts = {};
    int index = 0;
    int non_zero_counts = 0;
    //TODO: ask Jérôme if i should move that to constructor --> No
    for (std::vector<int> ix = gradients->new_index(); gradients->index_ok(ix); gradients->incr(ix)) {
      size_t count = gradients->samples->value(ix); 
      if (count > 0){
        insertIntoSortedList<size_t>(sorted_counts, count);
        non_zero_counts++;
      }
    }
    upper_threshold_count = sorted_counts[int(sorted_counts.size() * (1-lambda_max))];
    lower_threshold_count = sorted_counts[int(sorted_counts.size() * lambda_min)];
    // check this maybe we need to change before...
    sorted_counts.clear();
    int n_points = 0;
    for (std::vector<int> ix = gradients->new_index(); gradients->index_ok(ix); gradients->incr(ix)) {
      size_t count = gradients->samples->value(ix); 
      if (count < lower_threshold_count){
        sum_count += lower_threshold_count;
      }
      else if(count > upper_threshold_count){
        sum_count += upper_threshold_count;
      }
      else{
        sum_count += count;
      }
      n_points++;
    }
    m= float(sum_count)/n_points * 2/3;
    //TODO: ask Jérôme what is the difference between n_points and gradients->number_of_points() 
    std::cout << "m: " << m << " n_points: " << n_points << "gradients->number_of_points(): " << gradients->number_of_points() <<  std::endl;
    for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
        computation_grid->incr(ix)) {
      update_weighted_div_local(ix);
    }
    
  }

void integrate_potential::update_div_neighbors(const std::vector<int> &ix0)
{
  std::vector<int> ix(ix0);
  int i, j, k;

  // If not periodic, expanded grid ensures that upper neighbors of ix0 are valid grid points
  if (nd == 1) {
    return;

  } else if (nd == 2) {

    update_div_local(ix);
    ix[0]++; wrap(ix);
    update_div_local(ix);
    ix[1]++; wrap(ix);
    update_div_local(ix);
    ix[0]--; wrap(ix);
    update_div_local(ix);

  } else if (nd == 3) {

    for (i = 0; i<2; i++) {
      ix[1] = ix0[1];
      for (j = 0; j<2; j++) {
        ix[2] = ix0[2];
        for (k = 0; k<2; k++) {
          wrap(ix);
          update_div_local(ix);
          ix[2]++;
        }
        ix[1]++;
      }
      ix[0]++;
    }
  }
}



void integrate_potential::get_grad(cvm::real *g, std::vector<int> &ix)
{
  size_t i;
  bool edge = gradients->wrap_detect_edge(ix); // Detect edge if non-PBC

  if (edge) {
    for ( i = 0; i<nd; i++ ) {
          g[i] = 0.0;
    }
    return;
  }

  gradients->vector_value_smoothed(ix, g, b_smoothed);
}

void integrate_potential::update_div_local(const std::vector<int> &ix0)
{
  const size_t linear_index = address(ix0);
  int i, j, k;
  std::vector<int> ix = ix0;

  if (nd == 2) {
    // gradients at grid points surrounding the current scalar grid point
    cvm::real g00[2], g01[2], g10[2], g11[2];

    get_grad(g11, ix);
    ix[0] = ix0[0] - 1;
    get_grad(g01, ix);
    ix[1] = ix0[1] - 1;
    get_grad(g00, ix);
    ix[0] = ix0[0];
    get_grad(g10, ix);

    divergence[linear_index] = ((g10[0]-g00[0] + g11[0]-g01[0]) / widths[0]
                              + (g01[1]-g00[1] + g11[1]-g10[1]) / widths[1]) * 0.5;
  } else if (nd == 3) {
    cvm::real gc[24]; // stores 3d gradients in 8 contiguous bins
    int index = 0;

    ix[0] = ix0[0] - 1;
    for (i = 0; i<2; i++) {
      ix[1] = ix0[1] - 1;
      for (j = 0; j<2; j++) {
        ix[2] = ix0[2] - 1;
        for (k = 0; k<2; k++) {
          get_grad(gc + index, ix);
          index += 3;
          ix[2]++;
        }
        ix[1]++;
      }
      ix[0]++;
    }
  }
}

size_t integrate_potential::get_grad(std::vector<cvm::real> &g, std::vector<int> &ix){
  //TODO: it works fine
  size_t count = gradients->samples->value(ix);
  gradients -> vector_value(ix, g);
  return count;
}

inline size_t min(size_t a, size_t b) { return a < b ? a : b; }

void integrate_potential::prepare_divergence_calculation()
{
  surrounding_points_relative_positions.clear();
  int n_combinations = pow(2, nd);
  for (int i = 0; i < n_combinations; i++) {
    std::string binary = convert_base_two(i, nd);
    std::vector<int> surrounding_point_relative_position = {};
    for (char move_j : binary) {
      surrounding_point_relative_position.push_back(move_j - '0');
    }
    surrounding_points_relative_positions.push_back(surrounding_point_relative_position);
  }
}

void integrate_potential::update_weighted_div_local(const std::vector<int> &ix0)
/*
Updates the divergence at the point ix0
*/
{
  const size_t linear_index = computation_grid->address(ix0);
  int i, j, k;
  std::vector<int> ix = ix0;
  cvm::real div_at_point = 0;
  for (std::vector<int> surrounding_point_relative_position :
        surrounding_points_relative_positions) {
    std::vector<int> surrounding_point_coordinates = ix0;
    std::vector<cvm::real> gradient_at_surrounding_point(0, nd);
    for (int i = 0; i < nd; i++) {
      surrounding_point_coordinates[i] += surrounding_point_relative_position[i];
    }
    
    gradients->wrap_detect_edge(surrounding_point_coordinates);
    get_regularized_F(gradient_at_surrounding_point, surrounding_point_coordinates);
    cvm::real weight = get_regularized_weight(surrounding_point_coordinates);

    for (int i = 0; i < nd; i++) {
      div_at_point +=
           pow(-1, surrounding_point_relative_position[i] + 1) * gradient_at_surrounding_point[i] * weight / widths[i];
    }
  }
  
  divergence[linear_index] =
        div_at_point / pow(2, nd - 1);
  
}


/// Multiplication by sparse matrix representing Laplacian
/// NOTE: Laplacian must be symmetric for solving with CG
void integrate_potential::laplacian(const std::vector<cvm::real> &A, std::vector<cvm::real> &LA)
{
  if (nd == 2) {
    // DIMENSION 2

    size_t li, li2;
    int i, j;
    cvm::real fact;
    const cvm::real ffx = 1.0 / (widths[0] * widths[0]);
    const cvm::real ffy = 1.0 / (widths[1] * widths[1]);
    const int h = nx[1];
    const int w = nx[0];
    // offsets for 4 reference points of the Laplacian stencil
    int xm = -h;
    int xp = h;
    int ym = -1;
    int yp = 1;

    // NOTE on performance: this version is slightly sub-optimal because
    // it contains two double loops on the core of the array (for x and y terms)
    // The slightly faster version is in commit 0254cb5a2958cb2e135f268371c4b45fad34866b
    // yet it is much uglier, and probably horrible to extend to dimension 3
    // All terms in the matrix are assigned (=) during the x loops, then updated (+=)
    // with the y (and z) contributions


    // All x components except on x edges
    li = h; // Skip first column

    // Halve the term on y edges (if any) to preserve symmetry of the Laplacian matrix
    // (Long Chen, Finite Difference Methods, UCI, 2017)
    fact = periodic[1] ? 1.0 : 0.5;

    for (i = 1; i < w - 1; i++) {
      // Full range of j, but factor may change on y edges (j == 0 and j == h-1)
      LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
      li++;
      for (j = 1; j < h - 1; j++) {
        LA[li] = ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        li++;
      }
      LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
      li++;
    }
    // Edges along x (x components only)
    li = 0L;                              // Follows left edge
    li2 = h * static_cast<size_t>(w - 1); // Follows right edge
    if (periodic[0]) {
      xm = h * (w - 1);
      xp = h;
      fact = periodic[1] ? 1.0 : 0.5;
      LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
      LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
      li++;
      li2++;
      for (j = 1; j < h - 1; j++) {
        LA[li] = ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        LA[li2] = ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
        li++;
        li2++;
      }
      LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
      LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
    } else {
      xm = -h;
      xp = h;
      fact = periodic[1] ? 1.0 : 0.5; // Halve in corners in full PBC only
      // lower corner, "j == 0"
      LA[li] = fact * ffx * (A[li + xp] - A[li]);
      LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
      li++;
      li2++;
      for (j = 1; j < h - 1; j++) {
        // x gradient (+ y term of laplacian, calculated below)
        LA[li] = ffx * (A[li + xp] - A[li]);
        LA[li2] = ffx * (A[li2 + xm] - A[li2]);
        li++;
        li2++;
      }
      // upper corner, j == h-1
      LA[li] = fact * ffx * (A[li + xp] - A[li]);
      LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
    }

    // Now adding all y components
    // All y components except on y edges
    li = 1; // Skip first element (in first row)

    fact = periodic[0] ? 1.0 : 0.5; // for i == 0
    for (i = 0; i < w; i++) {
      // Factor of 1/2 on x edges if non-periodic
      if (i == 1)
        fact = 1.0;
      if (i == w - 1)
        fact = periodic[0] ? 1.0 : 0.5;
      for (j = 1; j < h - 1; j++) {
        LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        li++;
      }
      li += 2; // skip the edges and move to next column
    }
    // Edges along y (y components only)
    li = 0L;     // Follows bottom edge
    li2 = h - 1; // Follows top edge
    if (periodic[1]) {
      fact = periodic[0] ? 1.0 : 0.5;
      ym = h - 1;
      yp = 1;
      LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
      LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
      li += h;
      li2 += h;
      for (i = 1; i < w - 1; i++) {
        LA[li] += ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        LA[li2] += ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
        li += h;
        li2 += h;
      }
      LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
      LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
    } else {
      ym = -1;
      yp = 1;
      fact = periodic[0] ? 1.0 : 0.5; // Halve in corners in full PBC only
      // Left corner
      LA[li] += fact * ffy * (A[li + yp] - A[li]);
      LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
      li += h;
      li2 += h;
      for (i = 1; i < w - 1; i++) {
        // y gradient (+ x term of laplacian, calculated above)
        LA[li] += ffy * (A[li + yp] - A[li]);
        LA[li2] += ffy * (A[li2 + ym] - A[li2]);
        li += h;
        li2 += h;
      }
      // Right corner
      LA[li] += fact * ffy * (A[li + yp] - A[li]);
      LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
    }

  } else if (nd == 3) {
    // DIMENSION 3

    int i, j, k;
    size_t li, li2;
    cvm::real fact = 1.0;
    const cvm::real ffx = 1.0 / (widths[0] * widths[0]);
    const cvm::real ffy = 1.0 / (widths[1] * widths[1]);
    const cvm::real ffz = 1.0 / (widths[2] * widths[2]);
    const int h = nx[2]; // height
    const int d = nx[1]; // depth
    const int w = nx[0]; // width
    // offsets for 6 reference points of the Laplacian stencil
    int xm = -d * h;
    int xp = d * h;
    int ym = -h;
    int yp = h;
    int zm = -1;
    int zp = 1;

    cvm::real factx = periodic[0] ? 1 : 0.5; // factor to be applied on x edges
    cvm::real facty = periodic[1] ? 1 : 0.5; // same for y
    cvm::real factz = periodic[2] ? 1 : 0.5; // same for z
    cvm::real ifactx = 1 / factx;
    cvm::real ifacty = 1 / facty;
    cvm::real ifactz = 1 / factz;

    // All x components except on x edges
    li = d * static_cast<size_t>(h); // Skip left slab
    fact = facty * factz;
    for (i = 1; i < w - 1; i++) {
      for (j = 0; j < d; j++) { // full range of y
        if (j == 1)
          fact *= ifacty;
        if (j == d - 1)
          fact *= facty;
        LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        li++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) { // full range of z
          LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
          li++;
        }
        fact *= factz;
        LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        li++;
      }
    }
    // Edges along x (x components only)
    li = 0L;                                    // Follows left slab
    li2 = static_cast<size_t>(d) * h * (w - 1); // Follows right slab
    if (periodic[0]) {
      xm = d * h * (w - 1);
      xp = d * h;
      fact = facty * factz;
      for (j = 0; j < d; j++) {
        if (j == 1)
          fact *= ifacty;
        if (j == d - 1)
          fact *= facty;
        LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
        li++;
        li2++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
          LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
          li++;
          li2++;
        }
        fact *= factz;
        LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
        li++;
        li2++;
      }
    } else {
      xm = -d * h;
      xp = d * h;
      fact = facty * factz;
      for (j = 0; j < d; j++) {
        if (j == 1)
          fact *= ifacty;
        if (j == d - 1)
          fact *= facty;
        LA[li] = fact * ffx * (A[li + xp] - A[li]);
        LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
        li++;
        li2++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          // x gradient (+ y, z terms of laplacian, calculated below)
          LA[li] = fact * ffx * (A[li + xp] - A[li]);
          LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
          li++;
          li2++;
        }
        fact *= factz;
        LA[li] = fact * ffx * (A[li + xp] - A[li]);
        LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
        li++;
        li2++;
      }
    }

    // Now adding all y components
    // All y components except on y edges
    li = h; // Skip first column (in front slab)
    fact = factx * factz;
    for (i = 0; i < w; i++) { // full range of x
      if (i == 1)
        fact *= ifactx;
      if (i == w - 1)
        fact *= factx;
      for (j = 1; j < d - 1; j++) {
        LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        li++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
          li++;
        }
        fact *= factz;
        LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        li++;
      }
      li += 2 * h; // skip columns in front and back slabs
    }
    // Edges along y (y components only)
    li = 0L;                              // Follows front slab
    li2 = h * static_cast<size_t>(d - 1); // Follows back slab
    if (periodic[1]) {
      ym = h * (d - 1);
      yp = h;
      fact = factx * factz;
      for (i = 0; i < w; i++) {
        if (i == 1)
          fact *= ifactx;
        if (i == w - 1)
          fact *= factx;
        LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
        li++;
        li2++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
          LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
          li++;
          li2++;
        }
        fact *= factz;
        LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
        li++;
        li2++;
        li += h * static_cast<size_t>(d - 1);
        li2 += h * static_cast<size_t>(d - 1);
      }
    } else {
      ym = -h;
      yp = h;
      fact = factx * factz;
      for (i = 0; i < w; i++) {
        if (i == 1)
          fact *= ifactx;
        if (i == w - 1)
          fact *= factx;
        LA[li] += fact * ffy * (A[li + yp] - A[li]);
        LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
        li++;
        li2++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          // y gradient (+ x, z terms of laplacian, calculated above and below)
          LA[li] += fact * ffy * (A[li + yp] - A[li]);
          LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
          li++;
          li2++;
        }
        fact *= factz;
        LA[li] += fact * ffy * (A[li + yp] - A[li]);
        LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
        li++;
        li2++;
        li += h * static_cast<size_t>(d - 1);
        li2 += h * static_cast<size_t>(d - 1);
      }
    }

    // Now adding all z components
    // All z components except on z edges
    li = 1; // Skip first element (in bottom slab)
    fact = factx * facty;
    for (i = 0; i < w; i++) { // full range of x
      if (i == 1)
        fact *= ifactx;
      if (i == w - 1)
        fact *= factx;
      for (k = 1; k < h - 1; k++) {
        LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
        li++;
      }
      fact *= ifacty;
      li += 2;                      // skip edge slabs
      for (j = 1; j < d - 1; j++) { // full range of y
        for (k = 1; k < h - 1; k++) {
          LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
          li++;
        }
        li += 2; // skip edge slabs
      }
      fact *= facty;
      for (k = 1; k < h - 1; k++) {
        LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
        li++;
      }
      li += 2; // skip edge slabs
    }
    // Edges along z (z components onlz)
    li = 0;      // Follows bottom slab
    li2 = h - 1; // Follows top slab
    if (periodic[2]) {
      zm = h - 1;
      zp = 1;
      fact = factx * facty;
      for (i = 0; i < w; i++) {
        if (i == 1)
          fact *= ifactx;
        if (i == w - 1)
          fact *= factx;
        LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
        LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
        li += h;
        li2 += h;
        fact *= ifacty;
        for (j = 1; j < d - 1; j++) {
          LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
          LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
          li += h;
          li2 += h;
        }
        fact *= facty;
        LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
        LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
        li += h;
        li2 += h;
      }
    } else {
      zm = -1;
      zp = 1;
      fact = factx * facty;
      for (i = 0; i < w; i++) {
        if (i == 1)
          fact *= ifactx;
        if (i == w - 1)
          fact *= factx;
        LA[li] += fact * ffz * (A[li + zp] - A[li]);
        LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
        li += h;
        li2 += h;
        fact *= ifacty;
        for (j = 1; j < d - 1; j++) {
          // z gradient (+ x, y terms of laplacian, calculated above)
          LA[li] += fact * ffz * (A[li + zp] - A[li]);
          LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
          li += h;
          li2 += h;
        }
        fact *= facty;
        LA[li] += fact * ffz * (A[li + zp] - A[li]);
        LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
        li += h;
        li2 += h;
      }
    }
  }
}


/// Multiplication by sparse matrix representing Laplacian
/// NOTE: Laplacian must be symmetric for solving with CG
template<bool initialize_div_supplement> void integrate_potential::laplacian_weighted(const std::vector<cvm::real> &A, std::vector<cvm::real> &LA)
{
  for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix); computation_grid->incr(ix)){
    LA[computation_grid->address(ix)] = 0;
    if (initialize_div_supplement){
      div_border_supplement[computation_grid->address(ix)] = 0;
    }
  }
  // laplacian_matrix_test = std::vector<cvm::real>(computation_nt*computation_nt, 0);
  for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
        computation_grid->incr(ix)) {
      // TODO: delete this after testing
      // bool test = ix[0] == 128 && ix[1] == 118;
      for (std::pair<int, std::vector<int>> stencil_information: laplacian_stencil){
        std::vector<int> neighbor_relative_position = stencil_information.second;
        std::vector<int> neighbor_coordinate(nd, 0);
        for(int i = 0; i < nd; i++){
          neighbor_coordinate[i] = ix[i] + neighbor_relative_position[i];
        }
        bool virtual_point = computation_grid->wrap_detect_edge(neighbor_coordinate);
        cvm::real coefficient = calculate_weight_sum(neighbor_coordinate, weight_stencil[stencil_information.first])
                                * weight_counts[stencil_information.first] / pow(2, (nd-1)*2);
        std::pair<int, cvm::real> coefficient_regular_laplacian = neighbor_in_classic_laplacian_stencil[stencil_information.first];
        cvm::real coefficient_to_print = coefficient;

        
        if (coefficient_regular_laplacian.first){
          coefficient+= coefficient_regular_laplacian.second * m;
        }
        if (!virtual_point) {
          LA[computation_grid->address(ix)] += coefficient * A[computation_grid->address(neighbor_coordinate)];
          // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(neighbor_coordinate)] += coefficient;
        //   if (test){
        //   std::cout << "laplacian coordinates: " << "[" << computation_grid->address(ix) << ", " << computation_grid->address(neighbor_coordinate) << "]" << std::endl;
        //   std::cout << "coefficient: " << coefficient_to_print * pow(2, (nd-1)*2)<< std::endl;
        //   std::cout << "weight sum: " << calculate_weight_sum(neighbor_coordinate, weight_stencil[stencil_information.first]) << std::endl;
        //   std::cout << "weight counts: " << weight_counts[stencil_information.first] << std::endl;
        //   std::cout << "classical laplacian: " << coefficient_regular_laplacian.first << " " << coefficient_regular_laplacian.second << std::endl;
        //   std::cout << "virtual point: " << virtual_point << std::endl;
        //   std::cout << std::endl;
        // }
        } else {          
          std::vector<int> reference_point_coordinates(nd,0);
          computation_grid->wrap_to_edge(neighbor_coordinate, reference_point_coordinates);
          LA[computation_grid->address(ix)] += coefficient * A[computation_grid->address(reference_point_coordinates)];

          // if (test){
          //   std::cout << "laplacian coordinates: " << "[" << computation_grid->address(ix) << ", " << computation_grid->address(neighbor_coordinate) << "]" << std::endl;
          //   std::cout << "coefficient: " << coefficient_to_print * pow(2, (nd-1)*2)<< std::endl;
          //   std::cout << "weight sum: " << calculate_weight_sum(neighbor_coordinate, weight_stencil[stencil_information.first]) << std::endl;
          //   std::cout << "weight counts: " << weight_counts[stencil_information.first] << std::endl;
          //   std::cout << "classical laplacian: " << coefficient_regular_laplacian.first << " " << coefficient_regular_laplacian.second << std::endl;
          //   std::cout << "virtual point: " << virtual_point << std::endl;
          //   std::cout << "reference_point_coordinates: " << vec_to_string(reference_point_coordinates) << std::endl;
          //   std::cout << std::endl;
          // }
          // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(reference_point_coordinates)] += coefficient;

          cvm::real div_supplement_term = 0;
          if (initialize_div_supplement){
            std::vector<cvm::real> averaged_normal_vector = compute_averaged_border_normal_gradients(neighbor_coordinate);
            for (int i = 0; i < nd; i++){
              div_supplement_term += averaged_normal_vector[i] * neighbor_relative_position[i] * widths[i];
            }
          }
          bool test = reference_point_coordinates[0] == 128 && reference_point_coordinates[1] == 69;
          if (test){
            std::cout << "div_supplement_term: " << div_supplement_term << std::endl;
            std::cout << "coefficient: " << coefficient << std::endl;
          }
          div_border_supplement[computation_grid->address(ix)] -= div_supplement_term* coefficient;
        }
    }
  }
}

void integrate_potential::prepare_laplacian_calculation()
{
  for (int i = 0; i < std::pow(3, nd); i++) {
    std::string base_3 = convert_base_three(i);
    std::vector<int> direction;
    std::vector<std::vector<int>> weights_relative_positions = {
        {}}; // relative to the point of the stencil
    double weights_count = 0;
    int dim = 0;
    int number_of_non_zero_coordinates = 0;
    int non_zero_coordinate = -1;
    for (char direction_j : base_3) {
      int displacement_j = direction_j - '0';
      displacement_j -= 1;
      direction.push_back(displacement_j);
      switch (displacement_j) {
      case -1:
        weights_count += 1.0 / (widths[dim] * widths[dim]);
        weights_relative_positions =
            update_weight_relative_positions(weights_relative_positions, std::vector<int>{1});
        non_zero_coordinate = dim;
        number_of_non_zero_coordinates++;
        break;
      case 0:
        weights_count += -1.0 / (widths[dim] * widths[dim]);
        weights_relative_positions =
            update_weight_relative_positions(weights_relative_positions, std::vector<int>{0, 1});
        break;
      case 1:
        weights_count += 1.0 / (widths[dim] * widths[dim]);
        weights_relative_positions =
            update_weight_relative_positions(weights_relative_positions, std::vector<int>{0});
        non_zero_coordinate = dim;
        number_of_non_zero_coordinates++;
        break;
      }
      dim++;
    }
    // Store computed values in stencil maps
    laplacian_stencil[i] = direction;
    weight_stencil[i] = weights_relative_positions;
    weight_counts[i] = weights_count;

    // Store classic laplacian stencil information
    if (number_of_non_zero_coordinates <= 1) {
      if (non_zero_coordinate != -1)
        neighbor_in_classic_laplacian_stencil[i] = {
            true, 1 / (widths[non_zero_coordinate] * widths[non_zero_coordinate])};
      else {
        float sum = 0;
        for (int i = 0; i < nd; i++) {
          sum -= 2 / (widths[i] * widths[i]);
        }
        neighbor_in_classic_laplacian_stencil[i] = {true, sum};
      }
    } else {
      neighbor_in_classic_laplacian_stencil[i] = {false, -2};
    }
  }
}

void integrate_potential::print_laplacian_preparations()
{
  for (int i = 0; i < std::pow(3, nd); i++) {
    std::cout << "Stencil " << i << " is [";
    for (size_t j = 0; j < laplacian_stencil[i].size(); ++j) {
      std::cout << laplacian_stencil[i][j];
      if (j < laplacian_stencil[i].size() - 1)
        std::cout << ", ";
    }
    std::cout << "] with weights [";
  }
  std::cout << std::endl;
  std::cout << "weight stencil" << std::endl;
  for (int i = 0; i < std::pow(3, nd); i++) {
    std::cout << "Stencil " << i << " is [";
    for (size_t j = 0; j < weight_stencil[i].size(); ++j) {
      std::cout << vec_to_string(weight_stencil[i][j]);
      if (j < weight_stencil[i].size() - 1)
        std::cout << ", ";
    }
    std::cout << "]" << std::endl;
  }
  std::cout << std::endl;
  std::cout << "weight_counts" << std::endl;
  for (int i = 0; i < std::pow(3, nd); i++) {
    std::cout << "Stencil " << i << " is [";
    std::cout << weight_counts[i] << "[" << std::endl;
  }
  std::cout << std::endl;
  std::cout << "neighbor_in_classic_laplacian_stencil" << std::endl;
  for (int i = 0; i < std::pow(3, nd); i++) {
    std::cout << "Stencil " << i << " is [";
    std::cout << std::get<0>(neighbor_in_classic_laplacian_stencil[i]) << ", "
              << std::get<1>(neighbor_in_classic_laplacian_stencil[i]) << "]" << std::endl;
  }
  std::cout << std::endl;
}

std::vector<std::vector<int>> integrate_potential::update_weight_relative_positions(
    std::vector<std::vector<int>> &weights_relative_positions, std::vector<int> direction)
{

  std::vector<std::vector<int>> result;

  // For each weight direction and each existing relative position,
  // create a new position by appending the weight direction
  for (int weight_direction : direction) {
    for (const auto &weight_relative_position : weights_relative_positions) {
      // Clone the original position
      std::vector<int> weight_relative_position_clone = weight_relative_position;
      // Append the new direction
      weight_relative_position_clone.push_back(weight_direction);
      // Add to result
      result.push_back(weight_relative_position_clone);
    }
  }
  return result;
}


cvm::real integrate_potential::get_regularized_weight(std::vector<int> &ix){
  cvm::real regularized_weight;
  size_t count = gradients->samples->value(ix);
  if (count < lower_threshold_count)
  {
    regularized_weight = lower_threshold_count;
  }
  else if (count > upper_threshold_count){
    regularized_weight = upper_threshold_count;
  }
  else{
    regularized_weight = count;
  }
  return regularized_weight;

}

void integrate_potential::get_regularized_F(std::vector<cvm::real> &F, std::vector<int> &ix){
  // TODO: check if i cannot just use vector_value_smooth_instead/ old get_grad
  F.resize(nd);
  
  size_t count = get_grad(F, ix);
  float multiplier = 1;
  if (count < min_count_F){
    multiplier = 0;
  }
  else if (count < max_count_F){
    multiplier = (count - min_count_F) / (max_count_F - min_count_F); 
  }
  if (multiplier != 1){
    for (int i = 0; i < nd; i++){
      F[i] = multiplier * F[i];
    }
  }
}

cvm::real integrate_potential::calculate_weight_sum(std::vector<int> stencil_point,
                                                std::vector<std::vector<int>> directions)
{
  cvm::real weight_sum = 0;
  for (std::vector<int> direction : directions) {
    std::vector<int> weight_coordinate = stencil_point; // Initialize with stencil_point instead of size
    for (int i = 0; i < nd && i < direction.size(); i++) {
      weight_coordinate[i] += direction[i];
    }
    // bool test = stencil_point[0] == 128 && stencil_point[1] == 118;

    gradients->wrap_detect_edge(weight_coordinate);
    weight_sum += get_regularized_weight(weight_coordinate) - m;

    // if (test){
    //   std::cout << "weight_coordinate: " << vec_to_string(weight_coordinate) << std::endl;
    //   std::cout << "weight_sum: " << weight_sum << std::endl;
    // }
  }
  return weight_sum;
}

std::vector<cvm::real> integrate_potential::compute_averaged_border_normal_gradients(
    std::vector<int> virtual_point_coordinates)
{
  bool test = virtual_point_coordinates[0] == 129 && virtual_point_coordinates[1] == 69;
  std::vector<int> reference_point_coordinates(nd,0); // Initialize with correct size
  gradients->wrap_to_edge(virtual_point_coordinates, reference_point_coordinates);
  std::vector<int> directions_to_average_along;
  bool normal_directions[nd];
  for (int i = 0; i < nd; i++) {
    if ((0 <= virtual_point_coordinates[i] && virtual_point_coordinates[i] < computation_nx[i]) || periodic[i]) {
      directions_to_average_along.push_back(i);
      normal_directions[i] = false;
    } else {
      normal_directions[i] = true;
    }
  }
  // Find the positions of the gradients to average
  std::vector<std::vector<int>> gradients_to_average_relative_positions;
  if (directions_to_average_along.size() == 0) {
    std::vector<int> zero_vector(nd, 0);
    gradients_to_average_relative_positions.push_back(zero_vector);
  } else {
    for (int i = 0; i < pow(2, directions_to_average_along.size()); i++) {
      std::vector<int> gradient_to_average_relative_position(nd, 0);
      std::string binary = convert_base_two(i, directions_to_average_along.size());
      for (int bit_position = 0; bit_position < directions_to_average_along.size(); bit_position++) {
        gradient_to_average_relative_position[directions_to_average_along[bit_position]] =
            binary[bit_position] - '0';
      }
      gradients_to_average_relative_positions.push_back(gradient_to_average_relative_position);
    }
  }

  // compute the averaged bordered normal gradient
  std::vector<cvm::real> averaged_bordered_normal_gradient(nd, 0);
  // averaging the gradients
  for (int i = 0; i < gradients_to_average_relative_positions.size(); i++) {
    std::vector<int> gradient_position(reference_point_coordinates); // Initialize with reference_point_coordinates
    if (test){
      std::cout<< "gradient_position: " << vec_to_string(gradient_position) << std::endl;
      std::cout << vec_to_string(gradients_to_average_relative_positions[i]) << std::endl;
    }
    for (int j = 0; j < nd; j++) {
      gradient_position[j] += gradients_to_average_relative_positions[i][j];
    }
    std::vector<cvm::real> gradient(nd); // Initialize with correct size
    get_regularized_F(gradient, gradient_position);
    for (int j = 0; j < nd; j++) {
      averaged_bordered_normal_gradient[j] += gradient[j];
    }
  }
  // only keep the normal directions and average

  for (int j = 0; j < nd; j++) {
    if (!normal_directions[j]) {
      averaged_bordered_normal_gradient[j] = 0;
    }
    averaged_bordered_normal_gradient[j] /= gradients_to_average_relative_positions.size();
  }

  return averaged_bordered_normal_gradient;
}

std::string integrate_potential::convert_base_three(int n)
{
  std::string result = "";
  // Convert to base 3
  while (n > 0) {
    int remainder = n % 3;
    result.push_back('0' + remainder);
    n /= 3;
  }

  // Handle the case where n is 0
  if (result.empty()) {
    result = "0";
  }

  // Reverse the string (since we built it from right to left)
  std::reverse(result.begin(), result.end());

  // Pad with leading zeros if necessary
  while (result.size() < nd) {
    result = "0" + result;
  }

  // Truncate if the result has more digits than requested
  if (result.size() > nd) {
    result = result.substr(result.size() - nd);
  }
  return result;
}
std::string integrate_potential::convert_base_two(int n, int length)
{
  std::string result = "";

  // Convert to base 2
  while (n > 0) {
    int remainder = n % 2;
    result.push_back('0' + remainder);
    n /= 2;
  }

  // Handle the case where n is 0
  if (result.empty()) {
    result = "0";
  }

  // Reverse the string (since we built it from right to left)
  std::reverse(result.begin(), result.end());

  // Pad with leading zeros if necessary
  while (result.size() < length) {
    result = "0" + result;
  }

  // Truncate if the result has more digits than requested
  if (result.size() > length) {
    result = result.substr(result.size() - length);
  }
  return result;
}


void integrate_potential::nr_linbcg_sym(const bool weighted, const std::vector<cvm::real> &b,
                                        std::vector<cvm::real> &x, const cvm::real &tol,
                                        const int itmax, int &iter, cvm::real &err)
{
  cvm::real ak, akden, bk, bkden, bknum, bnrm;
  const cvm::real EPS = 1.0e-14;
  int j;
  std::vector<cvm::real> p(nt), r(nt), z(nt);
  typedef void (integrate_potential::*func_pointer)(const std::vector<double> &,
                                                      std::vector<double> &);
  func_pointer atimes =
      weighted ? &integrate_potential::laplacian_weighted<false> : &integrate_potential::laplacian;

  iter = 0;
  (this->*atimes)(x, r);
  for (j = 0; j < int(nt); j++) {
    r[j] = b[j] - r[j];
  }
  bnrm = l2norm(b);
  if (bnrm < EPS) {
    return; // Target is zero, will break relative error calc
  }
  //   asolve(r,z); // precon
  bkden = 1.0;
  while (iter < itmax) {
    ++iter;
    for (bknum = 0.0, j = 0; j < int(nt); j++) {
      bknum += r[j] * r[j]; // precon: z[j]*r[j]
    }
    if (iter == 1) {
      for (j = 0; j < int(nt); j++) {
        p[j] = r[j]; // precon: p[j] = z[j]
      }
    } else {
      bk = bknum / bkden;
      for (j = 0; j < int(nt); j++) {
        p[j] = bk * p[j] + r[j]; // precon:  bk*p[j] + z[j]
      }
    }
    bkden = bknum;
    (this->*atimes)(p, z);
    for (akden = 0.0, j = 0; j < int(nt); j++) {
      akden += z[j] * p[j];
    }
    ak = bknum / akden;
    for (j = 0; j < int(nt); j++) {
      x[j] += ak * p[j];
      r[j] -= ak * z[j];
    }
    //     asolve(r,z);  // precon
    err = l2norm(r) / bnrm;
    // if (cvm::debug())
      std::cout << "iter=" << std::setw(4) << iter + 1 << std::setw(12) << err << std::endl;
    if (err <= tol)
      break;
  }
}
template<typename T>
  typename std::vector<T>::iterator integrate_potential::insertIntoSortedList(std::vector<T>& sortedList, const T& value) {
    // Find the first position where the element is not less than value
    auto it = std::lower_bound(sortedList.begin(), sortedList.end(), value);
    
    // Insert the value at the found position and return iterator to inserted element
    return sortedList.insert(it, value);
}
cvm::real integrate_potential::l2norm(const std::vector<cvm::real> &x)
{
  size_t i;
  cvm::real sum = 0.0;
  for (i = 0; i < x.size(); i++)
    sum += x[i] * x[i];
  return sqrt(sum);
}
