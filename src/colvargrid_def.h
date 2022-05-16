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
  std::istream *is = cvm::main()->proxy->input_stream(filename, description);
  if (is == NULL) {
    return COLVARS_FILE_ERROR;
  }
  if (colvar_grid<T>::read_multicol(*is, add)) {
    cvm::main()->proxy->close_input_stream(filename);
    return COLVARS_OK;
  }
  return COLVARS_FILE_ERROR;
}


#endif
