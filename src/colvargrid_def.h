// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

/// \file Definition of the more complex members of colvar_grid<> template

#ifndef COLVARGRID_DEF_H
#define COLVARGRID_DEF_H

#include <iostream>
#include <iomanip>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvargrid.h"
#include "colvars_memstream.h"


template <class T, class IST> IST &read_restart_template(colvar_grid<T> &g, IST &is)
{
  auto const start_pos = is.tellg();
  std::string key, conf;
  if ((is >> key) && (key == std::string("grid_parameters"))) {
    is.seekg(start_pos);
    is >> colvarparse::read_block("grid_parameters", &conf);
    g.parse_params(conf, colvarparse::parse_silent);
  } else {
    cvm::log("Grid parameters are missing in the restart file, "
             "using those from the configuration.\n");
    is.seekg(start_pos);
  }
  g.read_raw(is);
  return is;
}


template <class T> std::istream &colvar_grid<T>::read_restart(std::istream &is)
{
  return read_restart_template<T, std::istream>(*this, is);
}


template <class T> cvm::memory_stream &colvar_grid<T>::read_restart(cvm::memory_stream &is)
{
  return read_restart_template<T, cvm::memory_stream>(*this, is);
}


template <class T> std::ostream &colvar_grid<T>::write_restart(std::ostream &os)
{
  write_params(os);
  write_raw(os);
  return os;
}


template <class T> cvm::memory_stream &colvar_grid<T>::write_restart(cvm::memory_stream &os)
{
  std::ostringstream oss;
  write_params(oss);
  os << oss.str();
  write_raw(os);
  return os;
}


template <class T, class IST> IST &read_raw_template(colvar_grid<T> &g, IST &is)
{
  auto const start_pos = is.tellg();

  for (std::vector<int> ix = g.new_index(); g.index_ok(ix); g.incr(ix)) {
    for (size_t imult = 0; imult < g.mult; imult++) {
      T new_value;
      if (is >> new_value) {
        g.value_input(ix, new_value, imult);
      } else {
        is.clear();
        is.seekg(start_pos);
        is.setstate(std::ios::failbit);
        cvm::error(
            "Error: failed to read all of the grid points from file.  Possible explanations: grid "
            "parameters in the configuration (lowerBoundary, upperBoundary, width) are different "
            "from those in the file, or the file is corrupt/incomplete.\n",
            COLVARS_INPUT_ERROR);
        return is;
      }
    }
  }

  g.has_data = true;
  return is;
}


template <class T> std::istream &colvar_grid<T>::read_raw(std::istream &is)
{
  return read_raw_template<T, std::istream>(*this, is);
}


template <class T> cvm::memory_stream &colvar_grid<T>::read_raw(cvm::memory_stream &is)
{
  return read_raw_template<T, cvm::memory_stream>(*this, is);
}


template <class T>
std::ostream &colvar_grid<T>::write_raw(std::ostream &os, size_t const buf_size) const
{
  auto const w = os.width();
  auto const p = os.precision();

  size_t count = 0;
  for (auto ix = new_index(); index_ok(ix); incr(ix)) {
    for (size_t imult = 0; imult < mult; imult++) {
      os << " " << std::setw(w) << std::setprecision(p) << value_output(ix, imult);
      if (((++count) % buf_size) == 0)
        os << "\n";
    }
  }
  // write a final newline only if buffer is not empty
  if ((count % buf_size) != 0)
    os << "\n";

  return os;
}


template <class T>
cvm::memory_stream &colvar_grid<T>::write_raw(cvm::memory_stream &os, size_t const buf_size) const
{
  for (auto ix = new_index(); index_ok(ix); incr(ix)) {
    for (size_t imult = 0; imult < mult; imult++) {
      os << value_output(ix, imult);
    }
  }
  return os;
}


template <class T>
std::istream & colvar_grid<T>::read_multicol(std::istream &is, bool add)
{
  // Data in the header: nColvars, then for each
  // xiMin, dXi, nPoints, periodic flag

  std::string   hash;
  cvm::real     lower, width, x;
  size_t        n, periodic_flag;
  bool          remap;
  std::vector<T>        new_value;
  std::vector<int>      nx_read;
  std::vector<int>      bin;

  if ( cv.size() > 0 && cv.size() != nd ) {
    cvm::error("Cannot read grid file: number of variables in file differs from number referenced by grid.\n");
    return is;
  }

  if ( !(is >> hash) || (hash != "#") ) {
    cvm::error("Error reading grid at position "+
               cvm::to_str(static_cast<size_t>(is.tellg()))+
               " in stream(read \"" + hash + "\")\n", COLVARS_INPUT_ERROR);
    return is;
  }

  is >> n;
  if ( n != nd ) {
    cvm::error("Error reading grid: wrong number of collective variables.\n");
    return is;
  }

  nx_read.resize(n);
  bin.resize(n);
  new_value.resize(mult);

  if (this->has_parent_data && add) {
    new_data.resize(data.size());
  }

  remap = false;
  for (size_t i = 0; i < nd; i++ ) {
    if ( !(is >> hash) || (hash != "#") ) {
      cvm::error("Error reading grid at position "+
                 cvm::to_str(static_cast<size_t>(is.tellg()))+
                 " in stream(read \"" + hash + "\")\n");
      return is;
    }

    is >> lower >> width >> nx_read[i] >> periodic_flag;


    if ( (cvm::fabs(lower - lower_boundaries[i].real_value) > 1.0e-10) ||
         (cvm::fabs(width - widths[i] ) > 1.0e-10) ||
         (nx_read[i] != nx[i]) ) {
      cvm::log("Warning: reading from different grid definition (colvar "
               + cvm::to_str(i+1) + "); remapping data on new grid.\n");
      remap = true;
    }
  }

  if ( remap ) {
    // re-grid data
    while (is.good()) {
      bool end_of_file = false;

      for (size_t i = 0; i < nd; i++ ) {
        if ( !(is >> x) ) end_of_file = true;
        bin[i] = value_to_bin_scalar(x, i);
      }
      if (end_of_file) break;

      for (size_t imult = 0; imult < mult; imult++) {
        is >> new_value[imult];
      }

      if ( index_ok(bin) ) {
        for (size_t imult = 0; imult < mult; imult++) {
          value_input(bin, new_value[imult], imult, add);
        }
      }
    }
  } else {
    // do not re-grid the data but assume the same grid is used
    for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix) ) {
      for (size_t i = 0; i < nd; i++ ) {
        is >> x;
      }
      for (size_t imult = 0; imult < mult; imult++) {
        is >> new_value[imult];
        value_input(ix, new_value[imult], imult, add);
      }
    }
  }
  has_data = true;
  return is;
}


template <class T>
int colvar_grid<T>::read_multicol(std::string const &filename,
                                  std::string description,
                                  bool add)
{
  std::istream &is = cvm::main()->proxy->input_stream(filename, description);
  if (!is) {
    return COLVARS_FILE_ERROR;
  }
  if (colvar_grid<T>::read_multicol(is, add)) {
    cvm::main()->proxy->close_input_stream(filename);
    return COLVARS_OK;
  }
  return COLVARS_FILE_ERROR;
}


template <class T>
std::ostream & colvar_grid<T>::write_multicol(std::ostream &os) const
{
  // Save the output formats
  std::ios_base::fmtflags prev_flags(os.flags());

  // Data in the header: nColvars, then for each
  // xiMin, dXi, nPoints, periodic

  os << std::setw(2) << "# " << nd << "\n";
  // Write the floating numbers in full precision
  os.setf(std::ios::scientific, std::ios::floatfield);
  for (size_t i = 0; i < nd; i++) {
    os << "# "
       << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << lower_boundaries[i] << " "
       << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << widths[i] << " "
       << std::setw(10) << nx[i] << "  "
       << periodic[i] << "\n";
  }

  for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix) ) {

    if (ix.back() == 0) {
      // if the last index is 0, add a new line to mark the new record
      os << "\n";
    }

    for (size_t i = 0; i < nd; i++) {
      os << " "
         << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec)
         << bin_to_value_scalar(ix[i], i);
    }
    os << " ";
    for (size_t imult = 0; imult < mult; imult++) {
      os << " "
         << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec)
         << value_output(ix, imult);
    }
    os << "\n";
  }

  // Restore the output formats
  os.flags(prev_flags);

  return os;
}


template <class T>
int colvar_grid<T>::write_multicol(std::string const &filename,
                                   std::string description) const
{
  int error_code = COLVARS_OK;
  std::ostream &os = cvm::main()->proxy->output_stream(filename, description);
  if (!os) {
    return COLVARS_FILE_ERROR;
  }
  error_code |= colvar_grid<T>::write_multicol(os) ? COLVARS_OK :
    COLVARS_FILE_ERROR;
  cvm::main()->proxy->close_output_stream(filename);
  return error_code;
}


template <class T>
std::ostream & colvar_grid<T>::write_opendx(std::ostream &os) const
{
  // write the header
  os << "object 1 class gridpositions counts";
  size_t icv;
  for (icv = 0; icv < num_variables(); icv++) {
    os << " " << number_of_points(icv);
  }
  os << "\n";

  os << "origin";
  for (icv = 0; icv < num_variables(); icv++) {
    os << " " << (lower_boundaries[icv].real_value + 0.5 * widths[icv]);
  }
  os << "\n";

  for (icv = 0; icv < num_variables(); icv++) {
    os << "delta";
    for (size_t icv2 = 0; icv2 < num_variables(); icv2++) {
      if (icv == icv2) os << " " << widths[icv];
      else os << " " << 0.0;
    }
    os << "\n";
  }

  os << "object 2 class gridconnections counts";
  for (icv = 0; icv < num_variables(); icv++) {
    os << " " << number_of_points(icv);
  }
  os << "\n";

  os << "object 3 class array type double rank 0 items "
     << number_of_points() << " data follows\n";

  write_raw(os);

  os << "object \"collective variables scalar field\" class field\n";
  return os;
}


template <class T>
int colvar_grid<T>::write_opendx(std::string const &filename,
                                 std::string description) const
{
  int error_code = COLVARS_OK;
  std::ostream &os = cvm::main()->proxy->output_stream(filename, description);
  if (!os) {
    return COLVARS_FILE_ERROR;
  }
  error_code |= colvar_grid<T>::write_opendx(os) ? COLVARS_OK :
    COLVARS_FILE_ERROR;
  cvm::main()->proxy->close_output_stream(filename);
  return error_code;
}

#endif
