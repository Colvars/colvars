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


template <class T>
colvar_grid<T>::colvar_grid(std::string const &filename, size_t mult_i)
{
  std::istream &is = cvm::main()->proxy->input_stream(filename, "multicol grid file");
  if (!is) {
    return;
  }

  std::string suffix = ".dx";
  if (filename.compare(filename.length() - suffix.length(),
                              suffix.length(),
                              suffix) == 0) {
    this->mult = mult_i;
    this->read_opendx(is);
  }
  else {
    std::string  hash;
    size_t i;

    if ( !(is >> hash) || (hash != "#") ) {
      cvm::error("Error reading grid at position "+
                  cvm::to_str(static_cast<size_t>(is.tellg()))+
                  " in stream(read \"" + hash + "\")\n");
      return;
    }

    is >> nd;
    mult = (mult_i == 0) ? nd : mult_i;

    std::vector<cvm::real> lower_in(nd), widths_in(nd);
    std::vector<int>       nx_in(nd);
    std::vector<int>       periodic_in(nd);

    for (i = 0; i < nd; i++ ) {
      if ( !(is >> hash) || (hash != "#") ) {
        cvm::error("Error reading grid at position "+
                    cvm::to_str(static_cast<size_t>(is.tellg()))+
                    " in stream(read \"" + hash + "\")\n");
        return;
      }

      is >> lower_in[i] >> widths_in[i] >> nx_in[i] >> periodic_in[i];
    }

    this->setup(nx_in, 0., mult);

    widths = widths_in;

    for (i = 0; i < nd; i++ ) {
      lower_boundaries.push_back(colvarvalue(lower_in[i]));
      periodic.push_back(static_cast<bool>(periodic_in[i]));
    }

    // Reset the istream for read_multicol, which expects the whole file
    is.clear();
    is.seekg(0);
    read_multicol(is);
  }
  cvm::main()->proxy->close_input_stream(filename);
  }


template <class T, class IST> IST &read_restart_template_(colvar_grid<T> &g, IST &is)
{
  auto const start_pos = is.tellg();
  std::string conf;
  if ((is >> colvarparse::read_block("grid_parameters", &conf)) &&
      (g.parse_params(conf, colvarparse::parse_restart) == COLVARS_OK) && g.read_raw(is)) {
    return is;
  }
  auto const error_pos = is.tellg();
  is.clear();
  is.seekg(start_pos);
  is.setstate(std::ios::failbit);
  cvm::error("Error: in reading grid state from stream at position " + cvm::to_str(error_pos) +
                 "\n",
             COLVARS_INPUT_ERROR);
  return is;
}


template <class T> std::istream &colvar_grid<T>::read_restart(std::istream &is)
{
  return read_restart_template_<T, std::istream>(*this, is);
}


template <class T> cvm::memory_stream &colvar_grid<T>::read_restart(cvm::memory_stream &is)
{
  return read_restart_template_<T, cvm::memory_stream>(*this, is);
}


template <class T> std::ostream &colvar_grid<T>::write_restart(std::ostream &os)
{
  os << "grid_parameters {\n" << get_state_params() << "}\n";
  write_raw(os);
  return os;
}


template <class T> cvm::memory_stream &colvar_grid<T>::write_restart(cvm::memory_stream &os)
{
  os << std::string("grid_parameters") << get_state_params();
  write_raw(os);
  return os;
}


template <class T, class IST> IST &read_raw_template_(colvar_grid<T> &g, IST &is)
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
  return read_raw_template_<T, std::istream>(*this, is);
}


template <class T> cvm::memory_stream &colvar_grid<T>::read_raw(cvm::memory_stream &is)
{
  return read_raw_template_<T, cvm::memory_stream>(*this, is);
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


template <class T> std::string colvar_grid<T>::get_state_params() const
{
  std::ostringstream os;
  size_t i;
  os << "  n_colvars " << nd << "\n";

  os << "  lower_boundaries ";
  for (i = 0; i < nd; i++)
    os << " " << lower_boundaries[i];
  os << "\n";

  os << "  upper_boundaries ";
  for (i = 0; i < nd; i++)
    os << " " << upper_boundaries[i];
  os << "\n";

  os << "  widths ";
  for (i = 0; i < nd; i++)
    os << " " << widths[i];
  os << "\n";

  os << "  sizes ";
  for (i = 0; i < nd; i++)
    os << " " << nx[i];
  os << "\n";

  return os.str();
}


template <class T> int colvar_grid<T>::parse_params(std::string const &conf,
                                                    colvarparse::Parse_Mode const parse_mode)
{
  if (cvm::debug())
    cvm::log("Reading grid configuration from string.\n");

  std::vector<int> old_nx = nx;
  std::vector<colvarvalue> old_lb = lower_boundaries;
  std::vector<colvarvalue> old_ub = upper_boundaries;
  std::vector<cvm::real> old_w = widths;

  {
    size_t nd_in = 0;
    // this is only used in state files
    colvarparse::get_keyval(conf, "n_colvars", nd_in, nd, colvarparse::parse_silent);
    if (nd_in != nd) {
      cvm::error("Error: trying to read data for a grid "
                 "that contains a different number of colvars ("+
                 cvm::to_str(nd_in)+") than the grid defined "
                 "in the configuration file("+cvm::to_str(nd)+
                 ").\n");
      return COLVARS_ERROR;
    }
  }

  // underscore keywords are used in state file
  colvarparse::get_keyval(conf, "lower_boundaries",
                          lower_boundaries, lower_boundaries, colvarparse::parse_silent);
  colvarparse::get_keyval(conf, "upper_boundaries",
                          upper_boundaries, upper_boundaries, colvarparse::parse_silent);
  // plural form is used in state file
  colvarparse::get_keyval(conf, "widths", widths, widths, colvarparse::parse_silent);

  // camel case keywords are used in config file
  colvarparse::get_keyval(conf, "lowerBoundary",
                          lower_boundaries, lower_boundaries, parse_mode);
  colvarparse::get_keyval(conf, "upperBoundary",
                          upper_boundaries, upper_boundaries, parse_mode);

  colvarparse::get_keyval(conf, "width", widths, widths, parse_mode);

  // only used in state file
  colvarparse::get_keyval(conf, "sizes", nx, nx, colvarparse::parse_silent);

  if (nd < lower_boundaries.size()) nd = lower_boundaries.size();

  if (! use_actual_value.size()) use_actual_value.assign(nd, false);
  if (! periodic.size()) periodic.assign(nd, false);
  if (! widths.size()) widths.assign(nd, 1.0);

  cvm::real eps = 1.e-10;

  bool new_params = false;
  if (old_nx.size()) {
    for (size_t i = 0; i < nd; i++) {
      if (old_nx[i] != nx[i] ||
          cvm::sqrt(cv[i]->dist2(old_lb[i], lower_boundaries[i])) > eps ||
          cvm::sqrt(cv[i]->dist2(old_ub[i], upper_boundaries[i])) > eps ||
          cvm::fabs(old_w[i] - widths[i]) > eps) {
        new_params = true;
      }
    }
  } else {
    new_params = true;
  }

  // reallocate the array in case the grid params have just changed
  if (new_params) {
    init_from_boundaries();
    // data.clear(); // no longer needed: setup calls clear()
    return this->setup(nx, T(), mult);
  }

  return COLVARS_OK;
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

  // Possible improvement:
  // - read line by line to ensure that layout is compatible with stated multiplicity
  // - detect end of file when not remapping

  if ( remap ) {
    // re-grid data
    while (is.good()) {
      bool end_of_file = false;

      for (size_t i = 0; i < nd; i++ ) {
        if ( !(is >> x) ) end_of_file = true;
        bin[i] = value_to_bin_scalar(x, i);
        // if x is out of bounds and we are using PBC, wrap it
        // Ignore out of bounds points in non-PBC
        wrap_detect_edge(bin);
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
std::istream & colvar_grid<T>::read_opendx(std::istream &is, bool add)
{
  std::string   line, token;
  std::vector<int> nx_read;
  std::vector<cvm::real> lower_in, widths_in;
  size_t        total_points = 1;
  size_t        num_items = 0;
  size_t        shape = 1;
  size_t        nd_read = 0;

  // 1. Parse OpenDX Header
  while (std::getline(is, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    iss >> token;

    if (token == "object") {
      iss >> token; // ID
      iss >> token; // "class"
      iss >> token; // class name

      if (token == "gridpositions") {
        iss >> token; // "counts"
        int n;
        while (iss >> n) {
          nx_read.push_back(n);
          total_points *= n;
        }
        nd_read = nx_read.size();
      }
      else if (token == "array") {
        while (iss >> token) {
          if (token == "shape") iss >> shape;
          if (token == "items") iss >> num_items;
          if (token == "data")  break;
        }
        num_items *= shape;
        break;
      }
    }
    else if (token == "origin") {
      cvm::real val;
      while (iss >> val) lower_in.push_back(val);
    }
    else if (token == "delta") {
      cvm::real val;
      std::vector<cvm::real> deltas;
      while (iss >> val) deltas.push_back(val);

      if (widths_in.size() < deltas.size()) {
        widths_in.push_back(deltas[widths_in.size()]);
      }
    }
  }

  if (nx_read.empty() || total_points == 0) {
    cvm::error("Error: invalid or empty OpenDX header.\n", COLVARS_INPUT_ERROR);
    return is;
  }

  size_t mult_read = num_items / total_points;
  if (mult_read == 0) mult_read = 1;

  bool remap = false;

  // 2. Setup the Grid OR Check Compatibility
  if (this->nx.empty()) {
    this->nd = nd_read;
    if (this->mult == 0) {
      this->mult = mult_read;
    }
    this->setup(nx_read, 0., this->mult);

    this->widths = widths_in;
    for (size_t i = 0; i < this->nd; i++) {
      // FIX: OpenDX origin is the center of the first bin.
      // colvar_grid lower_boundaries represents the lower edge of the first bin.
      // We subtract half the width to convert from node-centered to edge-centered.
      this->lower_boundaries.push_back(colvarvalue(lower_in[i] - 0.5 * widths_in[i]));
      this->periodic.push_back(false); // TODO: change that so it's not hard coded
    }
  } else {
    if (this->nd != nd_read) {
      cvm::error("Error reading OpenDX grid: wrong number of dimensions.\n");
      return is;
    }
    for (size_t i = 0; i < this->nd; i++) {
      // Reconstruct the expected OpenDX origin to accurately test compatibility
      cvm::real expected_dx_origin = this->lower_boundaries[i].real_value + 0.5 * this->widths[i];
      if ( (cvm::fabs(lower_in[i] - expected_dx_origin) > 1.0e-10) ||
           (cvm::fabs(widths_in[i] - this->widths[i]) > 1.0e-10) ||
           (nx_read[i] != this->nx[i]) ) {
        cvm::log("Warning: reading from different OpenDX grid definition; remapping data on new grid.\n");
        remap = true;
      }
    }
  }

  if (this->has_parent_data && add) {
    this->new_data.resize(this->data.size());
  }

  // 3. Read The Data Array
  std::vector<T> new_value(mult_read);
  std::vector<int> bin(this->nd);

  if (remap) {
    for (size_t pt = 0; pt < total_points; pt++) {
      std::vector<int> ix(this->nd, 0);
      size_t remainder = pt;

      for (int i = this->nd - 1; i >= 0; i--) {
        ix[i] = remainder % nx_read[i];
        remainder /= nx_read[i];
      }

      for (size_t imult = 0; imult < mult_read; imult++) {
        is >> new_value[imult];
      }

      for (size_t i = 0; i < this->nd; i++) {
        // lower_in holds the origin (bin center), mapping physically coordinates precisely
        cvm::real x = lower_in[i] + ix[i] * widths_in[i];
        bin[i] = this->value_to_bin_scalar(x, i);
      }
      this->wrap_detect_edge(bin);

      if (this->index_ok(bin)) {
        for (size_t imult = 0; imult < this->mult; imult++) {
          if (imult < mult_read) {
            this->value_input(bin, new_value[imult], imult, add);
          }
        }
      }
    }
  } else {
    for (std::vector<int> ix = this->new_index(); this->index_ok(ix); this->incr(ix)) {
      for (size_t imult = 0; imult < mult_read; imult++) {
        is >> new_value[imult];

        if (imult < this->mult) {
          this->value_input(ix, new_value[imult], imult, add);
        }
      }
    }
  }

  this->has_data = true;
  return is;
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
