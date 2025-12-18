// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARGRID_H
#define COLVARGRID_H

#include <algorithm>
#include <iosfwd>
#include <memory>

#include "colvar.h"
#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"

/// \brief Unified base class for grid of values of a function of several collective
/// variables
class colvar_grid_params  {

public:
  /// Number of dimensions
  size_t nd = 0;

  /// Number of points along each dimension
  std::vector<int> nx;

  /// Cumulative number of points along each dimension
  std::vector<int> nxc;

  /// Lower boundaries of the colvars in this grid
  std::vector<colvarvalue>  lower_boundaries;

  /// Upper boundaries of the colvars in this grid
  std::vector<colvarvalue>  upper_boundaries;

  /// Widths of the colvars in this grid
  std::vector<cvm::real>    widths;
};


/// \brief Grid of values of a function of several collective
/// variables \param T The data type
///
/// Only scalar colvars supported so far: vector colvars are treated as arrays
/// All common, type-independent members are collected in the base class colvar_grid_base
template <class T> class colvar_grid : public colvar_grid_params, public colvarparse {

  //protected:
public: // TODO create accessors for these after all instantiations work

  /// \brief Multiplicity of each datum (allow the binning of
  /// non-scalar types such as atomic gradients)
  size_t mult = 1;

  /// Total number of grid points
  size_t nt;

  /// Low-level array of values
  std::vector<T> data;

  /// Newly read data (used for count grids, when adding several grids read from disk)
  std::vector<T> new_data;

  /// Colvars collected in this grid
  std::vector<colvar *> cv;

  /// Do we request actual value (for extended-system colvars)?
  std::vector<bool> use_actual_value;

  /// Get the low-level (unrolled) address corresponding to an index
  inline size_t address(std::vector<int> const &ix) const
  {
    size_t addr = 0;
    for (size_t i = 0; i < nd; i++) {
      addr += ix[i]*static_cast<size_t>(nxc[i]);
      if (cvm::debug()) {
        if (ix[i] >= nx[i] || ix[i] < 0) {
          cvm::error("Error: index along dimension " + cvm::to_str(i) + " (" + cvm::to_str(ix[i])
            + ") exceeds bounds in colvar_grid.\n", COLVARS_BUG_ERROR);
          return 0;
        }
      }
    }
    return addr;
  }

  /// Inverse of address(): n-dim index from linear address
  inline void index(size_t address, std::vector<int> &ix) {
    for (int dim = nd - 1; dim >= 0; --dim) {
      ix[dim] = address % nx[dim];
      address /= nx[dim];
    }
  }

public:
  /// Whether some colvars are periodic
  std::vector<bool>        periodic;

  /// Whether some colvars have hard lower boundaries
  std::vector<bool>        hard_lower_boundaries;

  /// Whether some colvars have hard upper boundaries
  std::vector<bool>        hard_upper_boundaries;

  /// True if this is a count grid related to another grid of data
  bool has_parent_data;

  /// Whether this grid has been filled with data or is still empty
  bool has_data;

  /// Return the number of colvar objects
  inline size_t num_variables() const
  {
    return nd;
  }

  /// Return the numbers of points in all dimensions
  inline std::vector<int> const &number_of_points_vec() const
  {
    return nx;
  }

  /// Return the number of points in the i-th direction, if provided, or
  /// the total number
  inline size_t number_of_points(int const icv = -1) const
  {
    if (icv < 0) {
      return nt;
    } else {
      return nx[icv];
    }
  }

  /// Get the sizes in each direction
  inline std::vector<int> const & sizes() const
  {
    return nx;
  }

  /// Set the sizes in each direction
  inline void set_sizes(std::vector<int> const &new_sizes)
  {
    nx = new_sizes;
  }

  /// Return the multiplicity of the type used
  inline size_t multiplicity() const
  {
    return mult;
  }

  /// \brief Request grid to use actual values of extended coords
  inline void request_actual_value(bool b = true)
  {
    size_t i;
    for (i = 0; i < use_actual_value.size(); i++) {
      use_actual_value[i] = b;
    }
  }

  /// \brief Allocate data
  int setup(std::vector<int> const &nx_i,
            T const &t = T(),
            size_t const &mult_i = 1)
  {
    if (cvm::debug()) {
      cvm::log("Allocating grid: multiplicity = "+cvm::to_str(mult_i)+
               ", dimensions = "+cvm::to_str(nx_i)+".\n");
    }

    mult = mult_i;

    data.clear();

    nx = nx_i;
    nd = nx.size();

    nxc.resize(nd);

    // setup dimensions
    nt = mult;
    for (int i = nd-1; i >= 0; i--) {
      if (nx[i] <= 0) {
        cvm::error("Error: providing an invalid number of grid points, "+
                   cvm::to_str(nx[i])+".\n", COLVARS_BUG_ERROR);
        return COLVARS_ERROR;
      }
      nxc[i] = nt;
      nt *= nx[i];
    }

    if (cvm::debug()) {
      cvm::log("Total number of grid elements = "+cvm::to_str(nt)+".\n");
    }

    data.reserve(nt);
    data.assign(nt, t);

    return COLVARS_OK;
  }

  /// \brief Allocate data (allow initialization also after construction)
  int setup()
  {
    return setup(this->nx, T(), this->mult);
  }

  /// \brief Reset data (in case the grid is being reused)
  void reset(T const &t = T())
  {
    data.assign(nt, t);
  }


  /// Default constructor
  colvar_grid() : has_data(false)
  {
    nd = nt = 0;
    mult = 1;
    has_parent_data = false;
    this->setup();
  }

  /// Destructor
  virtual ~colvar_grid()
  {}

  /// \brief "Almost copy-constructor": only copies configuration
  /// parameters from another grid, but doesn't reallocate stuff;
  /// setup() must be called after that;
  colvar_grid(colvar_grid<T> const &g) : colvar_grid_params(colvar_grid_params(g)),
                                         colvarparse(),
                                         mult(g.mult),
                                         data(),
                                         cv(g.cv),
                                         use_actual_value(g.use_actual_value),
                                         periodic(g.periodic),
                                         hard_lower_boundaries(g.hard_lower_boundaries),
                                         hard_upper_boundaries(g.hard_upper_boundaries),
                                         has_parent_data(false),
                                         has_data(false)
  {}

  /// \brief Constructor from explicit grid sizes \param nx_i Number
  /// of grid points along each dimension \param t Initial value for
  /// the function at each point (optional) \param mult_i Multiplicity
  /// of each value
  colvar_grid(std::vector<int> const &nx_i,
              T const &t = T(),
              size_t mult_i = 1)
    : has_parent_data(false), has_data(false)
  {
    this->setup(nx_i, t, mult_i);
  }

  /// \brief Constructor from a vector of colvars or an optional grid config string
  /// \param add_extra_bin requests that non-periodic dimensions are extended
  /// by 1 bin to accommodate the integral (PMF) of another gridded quantity (gradient)
  colvar_grid(std::vector<colvar *> const &colvars,
              T const &t = T(),
              size_t mult_i = 1,
              bool add_extra_bin = false,
              std::shared_ptr<const colvar_grid_params> params = nullptr,
              std::string config = std::string())
    : has_parent_data(false), has_data(false)
  {
    (void) t;
    this->init_from_colvars(colvars, mult_i, add_extra_bin, params, config);
  }

  /// \brief Constructor from a multicol file
  /// \param filename multicol file containing data to be read
  /// \param multi_i multiplicity of the data - if 0, assume gradient multiplicity (mult = nd)
  colvar_grid(std::string const &filename, size_t mult_i = 1);

  int init_from_colvars(std::vector<colvar *> const &colvars,
                        size_t mult_i = 1,
                        bool add_extra_bin = false,
                        std::shared_ptr<const colvar_grid_params> params = nullptr,
                        std::string config = std::string())
  {
    if (cvm::debug()) {
      cvm::log("Reading grid configuration from collective variables.\n");
    }

    cv = colvars;
    nd = colvars.size();
    mult = mult_i;

    size_t i;

    if (cvm::debug()) {
      cvm::log("Allocating a grid for "+cvm::to_str(colvars.size())+
               " collective variables, multiplicity = "+cvm::to_str(mult_i)+".\n");
    }

    for (i =  0; i < nd; i++) {
      if (cv[i]->value().type() != colvarvalue::type_scalar) {
        cvm::error("Colvar grids can only be automatically "
                   "constructed for scalar variables.  "
                   "ABF and histogram can not be used; metadynamics "
                   "can be used with useGrids disabled.\n", COLVARS_INPUT_ERROR);
        return COLVARS_ERROR;
      }

      if (cv[i]->width <= 0.0) {
        cvm::error("Tried to initialize a grid on a "
                   "variable with negative or zero width.\n", COLVARS_INPUT_ERROR);
        return COLVARS_ERROR;
      }

      widths.push_back(cv[i]->width);
      hard_lower_boundaries.push_back(cv[i]->is_enabled(colvardeps::f_cv_hard_lower_boundary));
      hard_upper_boundaries.push_back(cv[i]->is_enabled(colvardeps::f_cv_hard_upper_boundary));

      // By default, get reported colvar value (for extended Lagrangian colvars)
      use_actual_value.push_back(false);

      // except if a colvar is specified twice in a row
      // then the first instance is the actual value
      // For histograms of extended-system coordinates
      if (i > 0 && cv[i-1] == cv[i]) {
        use_actual_value[i-1] = true;
      }

      // This needs to work if the boundaries are undefined in the colvars
      lower_boundaries.push_back(cv[i]->lower_boundary);
      upper_boundaries.push_back(cv[i]->upper_boundary);
    }

    // Replace widths and boundaries with optional custom configuration
    if (!config.empty()) {
      this->parse_params(config);
      this->check_keywords(config, "grid");

      if (params) {
        cvm::error("Error: init_from_colvars was passed both a grid config and a template grid.", COLVARS_BUG_ERROR);
        return COLVARS_BUG_ERROR;
      }
    } else if (params) {
      // Match grid sizes with template

      if (params->nd != nd) {
        cvm::error("Trying to initialize grid from template with wrong dimension (" +
                    cvm::to_str(params->nd) + " instead of " +
                    cvm::to_str(this->nd) + ").");
        return COLVARS_ERROR;
      }

      widths =params->widths;
      lower_boundaries =params->lower_boundaries;
      upper_boundaries =params->upper_boundaries;
    }

    // Only now can we determine periodicity
    for (i =  0; i < nd; i++) {
      periodic.push_back(cv[i]->periodic_boundaries(lower_boundaries[i].real_value,
                                                    upper_boundaries[i].real_value));

      if (add_extra_bin) {
        // Shift the grid by half the bin width (values at edges instead of center of bins)
        lower_boundaries[i] -= 0.5 * widths[i];

        if (periodic[i]) {
          // Just shift
          upper_boundaries[i] -= 0.5 * widths[i];
        } else {
          // Widen grid by one bin width
          upper_boundaries[i] += 0.5 * widths[i];
        }
      }
    }

    // Reset grid sizes based on widths and boundaries
    this->init_from_boundaries();
    return this->setup();
  }

  int init_from_boundaries()
  {
    if (cvm::debug()) {
      cvm::log("Configuring grid dimensions from colvars boundaries.\n");
    }

    // these will have to be updated
    nx.clear();
    nxc.clear();
    nt = 0;

    for (size_t i =  0; i < lower_boundaries.size(); i++) {
      // Re-compute periodicity using current grid boundaries
      periodic[i] = cv[i]->periodic_boundaries(lower_boundaries[i].real_value,
                                               upper_boundaries[i].real_value);

      cvm::real nbins = ( upper_boundaries[i].real_value -
                          lower_boundaries[i].real_value ) / widths[i];
      int nbins_round = (int)(nbins+0.5);

      if (cvm::fabs(nbins - cvm::real(nbins_round)) > 1.0E-10) {
        cvm::log("Warning: grid interval("+
                 cvm::to_str(lower_boundaries[i], cvm::cv_width, cvm::cv_prec)+" - "+
                 cvm::to_str(upper_boundaries[i], cvm::cv_width, cvm::cv_prec)+
                 ") is not commensurate to its bin width("+
                 cvm::to_str(widths[i], cvm::cv_width, cvm::cv_prec)+").\n");
        upper_boundaries[i].real_value = lower_boundaries[i].real_value +
          (nbins_round * widths[i]);
      }

      if (cvm::debug())
        cvm::log("Number of points is "+cvm::to_str((int) nbins_round)+
                 " for the colvar no. "+cvm::to_str(i+1)+".\n");

      nx.push_back(nbins_round);
    }

    return COLVARS_OK;
  }

  /// Wrap an index vector around periodic boundary conditions
  /// also checks validity of non-periodic indices
  inline void wrap(std::vector<int> & ix) const
  {
    for (size_t i = 0; i < nd; i++) {
      if (periodic[i]) {
        ix[i] = (ix[i] + nx[i]) % nx[i]; // Avoid modulo with negative operands (implementation-defined)
      } else {
        if (ix[i] < 0 || ix[i] >= nx[i]) {
          cvm::error("Trying to wrap illegal index vector (non-PBC) for a grid point: "
                     + cvm::to_str(ix), COLVARS_BUG_ERROR);
          return;
        }
      }
    }
  }

  /// Wrap an index vector around periodic boundary conditions
  /// or detects edges if non-periodic
  inline bool wrap_detect_edge(std::vector<int> & ix) const
  {
    bool edge = false;
    for (size_t i = 0; i < nd; i++) {
      if (periodic[i]) {
        ix[i] = (ix[i] + nx[i]) % nx[i]; // Avoid modulo with negative operands (implementation-defined)
      } else if (ix[i] < 0 || ix[i] >= nx[i]) {
        edge = true;
      }
    }
    return edge;
  }

  /// Wrap an index vector around periodic boundary conditions
  /// or brings back to nearest edge if non-periodic
  inline bool wrap_to_edge(std::vector<int> & ix, std::vector<int> & edge_bin) const
  {
    bool edge = false;
    edge_bin = ix;
    for (size_t i = 0; i < nd; i++) {
      if (periodic[i]) {
        ix[i] = (ix[i] + nx[i]) % nx[i]; // Avoid modulo with negative operands (implementation-defined)
        edge_bin[i] = ix[i];
      } else if (ix[i] < 0) {
        edge = true;
        edge_bin[i] = 0;
      } else if (ix[i] >= nx[i]) {
        edge = true;
        edge_bin[i] = nx[i] - 1;
      }
    }
    return edge;
  }


  /// \brief Report the bin corresponding to the current value of variable i
  inline int current_bin_scalar(int const i) const
  {
    return value_to_bin_scalar(use_actual_value[i] ? cv[i]->actual_value() : cv[i]->value(), i);
  }

  /// \brief Report the flattened bin address corresponding to the current value of all variables
  /// and assign first or last bin if out of boundaries
  inline int current_bin_flat_bound() const
  {
    std::vector<int> index = new_index();
    for (size_t i = 0; i < nd; i++) {
      index[i] = current_bin_scalar_bound(i);
    }
    return address(index);
  }

  /// \brief Report the bin corresponding to the current value of variable i
  /// and assign first or last bin if out of boundaries
  inline int current_bin_scalar_bound(int const i) const
  {
    return value_to_bin_scalar_bound(use_actual_value[i] ? cv[i]->actual_value() : cv[i]->value(), i);
  }

  /// \brief Report the bin corresponding to the current value of item iv in variable i
  inline int current_bin_scalar(int const i, int const iv) const
  {
    return value_to_bin_scalar(use_actual_value[i] ?
                               cv[i]->actual_value().vector1d_value[iv] :
                               cv[i]->value().vector1d_value[iv], i);
  }

  /// \brief Use the lower boundary and the width to report which bin
  /// the provided value is in
  inline int value_to_bin_scalar(colvarvalue const &value, const int i) const
  {
    return (int) cvm::floor( (value.real_value - lower_boundaries[i].real_value) / widths[i] );
  }

  /// \brief Report the fraction of bin beyond current_bin_scalar()
  inline cvm::real current_bin_scalar_fraction(int const i) const
  {
    return value_to_bin_scalar_fraction(use_actual_value[i] ? cv[i]->actual_value() : cv[i]->value(), i);
  }

  /// \brief Use the lower boundary and the width to report the fraction of bin
  /// beyond value_to_bin_scalar() that the provided value is in
  inline cvm::real value_to_bin_scalar_fraction(colvarvalue const &value, const int i) const
  {
    cvm::real x = (value.real_value - lower_boundaries[i].real_value) / widths[i];
    return x - cvm::floor(x);
  }

  /// \brief Use the lower boundary and the width to report which bin
  /// the provided value is in and assign first or last bin if out of boundaries
  inline int value_to_bin_scalar_bound(colvarvalue const &value, const int i) const
  {
    int bin_index = cvm::floor( (value.real_value - lower_boundaries[i].real_value) / widths[i] );

    // Wrap bins for periodic dimensions before truncating
    if (periodic[i]) bin_index %= nx[i];
    if (bin_index < 0) bin_index=0;
    if (bin_index >=int(nx[i])) bin_index=int(nx[i])-1;
    return (int) bin_index;
  }

  /// \brief Same as the standard version, but uses another grid definition
  inline int value_to_bin_scalar(colvarvalue const &value,
                                 colvarvalue const &new_offset,
                                 cvm::real   const &new_width) const
  {
    return (int) cvm::floor( (value.real_value - new_offset.real_value) / new_width );
  }

  /// \brief Use the two boundaries and the width to report the
  /// central value corresponding to a bin index
  inline colvarvalue bin_to_value_scalar(int const &i_bin, int const i) const
  {
    return lower_boundaries[i].real_value + widths[i] * (0.5 + i_bin);
  }

  /// \brief Same as the standard version, but uses different parameters
  inline colvarvalue bin_to_value_scalar(int const &i_bin,
                                         colvarvalue const &new_offset,
                                         cvm::real const &new_width) const
  {
    return new_offset.real_value + new_width * (0.5 + i_bin);
  }

  /// Set the value at the point with index ix
  inline void set_value(std::vector<int> const &ix,
                        T const &t,
                        size_t const &imult = 0)
  {
    data[this->address(ix)+imult] = t;
    has_data = true;
  }

  /// Set the value at the point with linear address i (for speed)
  inline void set_value(size_t i, T const &t)
  {
    data[i] = t;
  }

 /// Get the value at the point with linear address i (for speed)
  inline T get_value(size_t i) const
  {
    return data[i];
  }


  /// \brief Get the change from this to other_grid
  /// and store the result in this.
  /// this_grid := other_grid - this_grid
  /// Grids must have the same dimensions.
  void delta_grid(colvar_grid<T> const &other_grid)
  {

    if (other_grid.multiplicity() != this->multiplicity()) {
      cvm::error("Error: trying to subtract two grids with "
                 "different multiplicity.\n");
      return;
    }

    if (other_grid.data.size() != this->data.size()) {
      cvm::error("Error: trying to subtract two grids with "
                 "different size.\n");
      return;
    }

    for (size_t i = 0; i < data.size(); i++) {
      data[i] = other_grid.data[i] - data[i];
    }
    has_data = true;
  }


  /// \brief Copy data from another grid of the same type, AND
  /// identical definition (boundaries, widths)
  /// Added for shared ABF.
  void copy_grid(colvar_grid<T> const &other_grid)
  {
    if (other_grid.multiplicity() != this->multiplicity()) {
      cvm::error("Error: trying to copy two grids with "
                 "different multiplicity.\n");
      return;
    }

    if (other_grid.data.size() != this->data.size()) {
      cvm::error("Error: trying to copy two grids with "
                 "different size.\n");
      return;
    }


    for (size_t i = 0; i < data.size(); i++) {
      data[i] = other_grid.data[i];
    }
    has_data = true;
  }

  /// \brief Extract the grid data as they are represented in memory.
  /// Put the results in "out_data".
  void raw_data_out(T* out_data) const
  {
    for (size_t i = 0; i < data.size(); i++) out_data[i] = data[i];
  }
  void raw_data_out(std::vector<T>& out_data) const
  {
    out_data = data;
  }
  /// \brief Input the data as they are represented in memory.
  void raw_data_in(const T* in_data)
  {
    for (size_t i = 0; i < data.size(); i++) data[i] = in_data[i];
    has_data = true;
  }
  void raw_data_in(const std::vector<T>& in_data)
  {
    data = in_data;
    has_data = true;
  }
  /// \brief Size of the data as they are represented in memory.
  size_t raw_data_num() const { return data.size(); }


  /// \brief Get the binned value indexed by ix, or the first of them
  /// if the multiplicity is larger than 1
  inline T const & value(std::vector<int> const &ix,
                         size_t const &imult = 0) const
  {
    return data[this->address(ix) + imult];
  }

  /// \brief Get the binned value indexed by linear address i
  inline T const & value(size_t i) const
  {
    return data[i];
  }

  /// \brief Add a constant to all elements (fast loop)
  inline void add_constant(T const &t)
  {
    for (size_t i = 0; i < nt; i++)
      data[i] += t;
    has_data = true;
  }

  /// \brief Multiply all elements by a scalar constant (fast loop)
  inline void multiply_constant(cvm::real const &a)
  {
    for (size_t i = 0; i < nt; i++)
      data[i] *= a;
  }

  /// \brief Assign values that are smaller than scalar constant the latter value (fast loop)
  inline void remove_small_values(cvm::real const &a)
  {
    for (size_t i = 0; i < nt; i++)
      if(data[i]<a) data[i] = a;
  }


  /// \brief Get the bin indices corresponding to the provided values of
  /// the colvars
  inline std::vector<int> const get_colvars_index(std::vector<colvarvalue> const &values) const
  {
    std::vector<int> index = new_index();
    for (size_t i = 0; i < nd; i++) {
      index[i] = value_to_bin_scalar(values[i], i);
    }
    return index;
  }

  /// \brief Get the bin indices corresponding to the current values
  /// of the colvars
  inline std::vector<int> const get_colvars_index() const
  {
    std::vector<int> index = new_index();
    for (size_t i = 0; i < nd; i++) {
      index[i] = current_bin_scalar(i);
    }
    return index;
  }

  /// \brief Get the bin indices corresponding to the provided values of
  /// the colvars and assign first or last bin if out of boundaries
  inline std::vector<int> const get_colvars_index_bound() const
  {
    std::vector<int> index = new_index();
    for (size_t i = 0; i < nd; i++) {
      index[i] = current_bin_scalar_bound(i);
    }
    return index;
  }

  /// \brief Get the minimal distance (in number of bins) from the
  /// boundaries; a negative number is returned if the given point is
  /// off-grid
  inline cvm::real bin_distance_from_boundaries(std::vector<colvarvalue> const &values,
                                                bool skip_hard_boundaries = false)
  {
    cvm::real minimum = 1.0E+16;
    for (size_t i = 0; i < nd; i++) {

      if (periodic[i]) continue;

      cvm::real dl = cvm::sqrt(cv[i]->dist2(values[i], lower_boundaries[i])) / widths[i];
      cvm::real du = cvm::sqrt(cv[i]->dist2(values[i], upper_boundaries[i])) / widths[i];

      if (values[i].real_value < lower_boundaries[i])
        dl *= -1.0;
      if (values[i].real_value > upper_boundaries[i])
        du *= -1.0;

      if ( ((!skip_hard_boundaries) || (!hard_lower_boundaries[i])) && (dl < minimum))
        minimum = dl;
      if ( ((!skip_hard_boundaries) || (!hard_upper_boundaries[i])) && (du < minimum))
        minimum = du;
    }

    return minimum;
  }


  /// \brief Add data from another grid of the same type
  ///
  /// Note: this function maps other_grid inside this one regardless
  /// of whether it fits or not.
  void map_grid(colvar_grid<T> const &other_grid)
  {
    if (other_grid.multiplicity() != this->multiplicity()) {
      cvm::error("Error: trying to merge two grids with values of "
                 "different multiplicity.\n");
      return;
    }

    std::vector<colvarvalue> const &gb  = this->lower_boundaries;
    std::vector<cvm::real> const &gw    = this->widths;
    std::vector<colvarvalue> const &ogb = other_grid.lower_boundaries;
    std::vector<cvm::real> const &ogw   = other_grid.widths;

    std::vector<int> ix = this->new_index();
    std::vector<int> oix = other_grid.new_index();

    if (cvm::debug())
      cvm::log("Remapping grid...\n");
    for ( ; this->index_ok(ix); this->incr(ix)) {

      for (size_t i = 0; i < nd; i++) {
        oix[i] =
          value_to_bin_scalar(bin_to_value_scalar(ix[i], gb[i], gw[i]),
                              ogb[i],
                              ogw[i]);
      }

      if (! other_grid.index_ok(oix)) {
        continue;
      }

      for (size_t im = 0; im < mult; im++) {
        this->set_value(ix, other_grid.value(oix, im), im);
      }
    }

    has_data = true;
    if (cvm::debug())
      cvm::log("Remapping done.\n");
  }

  /// \brief Add data from another grid of the same type, AND
  /// identical definition (boundaries, widths)
  void add_grid(colvar_grid<T> const &other_grid,
                cvm::real scale_factor = 1.0)
  {
    if (other_grid.multiplicity() != this->multiplicity()) {
      cvm::error("Error: trying to sum togetehr two grids with values of "
                 "different multiplicity.\n");
      return;
    }
    if (scale_factor != 1.0)
      for (size_t i = 0; i < data.size(); i++) {
        data[i] += static_cast<T>(scale_factor * other_grid.data[i]);
      }
    else
      // skip multiplication if possible
      for (size_t i = 0; i < data.size(); i++) {
        data[i] += other_grid.data[i];
      }
    has_data = true;
  }

  /// \brief Return the value suitable for output purposes (so that it
  /// may be rescaled or manipulated without changing it permanently)
  virtual T value_output(std::vector<int> const &ix,
                         size_t const &imult = 0) const
  {
    return value(ix, imult);
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (the two may be different,
  /// e.g. when using colvar_grid_count)
  virtual void value_input(std::vector<int> const &ix,
                           T const &t,
                           size_t const &imult = 0,
                           bool add = false)
  {
    if ( add )
      data[address(ix) + imult] += t;
    else
      data[address(ix) + imult] = t;
    has_data = true;
  }


  //   /// Get the pointer to the binned value indexed by ix
  //   inline T const *value_p (std::vector<int> const &ix)
  //   {
  //     return &(data[address (ix)]);
  //   }

  /// \brief Get the index corresponding to the "first" bin, to be
  /// used as the initial value for an index in looping
  inline std::vector<int> const new_index() const
  {
    return std::vector<int> (nd, 0);
  }

  /// \brief Check that the index is within range in each of the
  /// dimensions
  inline bool index_ok(std::vector<int> const &ix) const
  {
    for (size_t i = 0; i < nd; i++) {
      if ( (ix[i] < 0) || (ix[i] >= int(nx[i])) )
        return false;
    }
    return true;
  }

  /// \brief Increment the index, in a way that will make it loop over
  /// the whole nd-dimensional array
  inline void incr(std::vector<int> &ix) const
  {
    for (int i = ix.size()-1; i >= 0; i--) {

      ix[i]++;

      if (ix[i] >= nx[i]) {

        if (i > 0) {
          ix[i] = 0;
          continue;
        } else {
          // this is the last iteration, a non-valid index is being
          // set for the outer index, which will be caught by
          // index_ok()
          ix[0] = nx[0];
          return;
        }
      } else {
        return;
      }
    }
  }

  /// Write the current grid parameters to a string
  virtual std::string get_state_params() const;

  /// Read new grid parameters from a string
  int parse_params(std::string const &conf,
                   colvarparse::Parse_Mode const parse_mode = colvarparse::parse_normal);

  /// \brief Check that the grid information inside (boundaries,
  /// widths, ...) is consistent with the current setting of the
  /// colvars
  void check_consistency()
  {
    for (size_t i = 0; i < nd; i++) {
      if ( (cvm::sqrt(cv[i]->dist2(cv[i]->lower_boundary,
                                   lower_boundaries[i])) > 1.0E-10) ||
           (cvm::sqrt(cv[i]->dist2(cv[i]->upper_boundary,
                                   upper_boundaries[i])) > 1.0E-10) ||
           (cvm::sqrt(cv[i]->dist2(cv[i]->width,
                                   widths[i])) > 1.0E-10) ) {
        cvm::error("Error: restart information for a grid is "
                   "inconsistent with that of its colvars.\n");
        return;
      }
    }
  }


  /// \brief Check that the grid information inside (boundaries,
  /// widths, ...) is consistent with that of another grid
  void check_consistency(colvar_grid<T> const &other_grid)
  {
    for (size_t i = 0; i < nd; i++) {
      // we skip dist2(), because periodicities and the like should
      // matter: boundaries should be EXACTLY the same (otherwise,
      // map_grid() should be used)
      if ( (cvm::fabs(other_grid.lower_boundaries[i] -
                      lower_boundaries[i]) > 1.0E-10) ||
           (cvm::fabs(other_grid.upper_boundaries[i] -
                      upper_boundaries[i]) > 1.0E-10) ||
           (cvm::fabs(other_grid.widths[i] -
                      widths[i]) > 1.0E-10) ||
           (data.size() != other_grid.data.size()) ) {
        cvm::error("Error: inconsistency between "
                   "two grids that are supposed to be equal, "
                   "aside from the data stored.\n");
        return;
      }
    }
  }

  /// Read all grid parameters and data from a formatted stream
  std::istream & read_restart(std::istream &is);

  /// Read all grid parameters and data from an unformatted stream
  cvm::memory_stream & read_restart(cvm::memory_stream &is);

  /// Write all grid parameters and data to a formatted stream
  std::ostream & write_restart(std::ostream &os);

  /// Write all grid parameters and data to an unformatted stream
  cvm::memory_stream & write_restart(cvm::memory_stream &os);

  /// Read all grid parameters and data from a formatted stream
  std::istream &read_raw(std::istream &is);

  /// Read all grid parameters and data from an unformatted stream
  cvm::memory_stream &read_raw(cvm::memory_stream &is);

  /// Write all grid data to a formatted stream (without labels, as they are represented in memory)
  /// \param[in,out] os Stream object
  /// \param[in] buf_size Number of values per line
  std::ostream &write_raw(std::ostream &os, size_t const buf_size = 3) const;

  /// Write all grid data to an unformatted stream
  /// \param[in,out] os Stream object
  /// \param[in] buf_size Number of values per line (note: ignored because there is no formatting)
  cvm::memory_stream &write_raw(cvm::memory_stream &os, size_t const buf_size = 3) const;

  /// Read a grid written by write_multicol(), incrementing if add is true
  std::istream & read_multicol(std::istream &is, bool add = false);

  /// Read a grid written by write_multicol(), incrementing if add is true
  int read_multicol(std::string const &filename,
                    std::string description = "grid file",
                    bool add = false);

  /// Write grid in a format which is both human-readable and gnuplot-friendly
  std::ostream & write_multicol(std::ostream &os) const;

  /// Write grid in a format which is both human-readable and gnuplot-friendly
  int write_multicol(std::string const &filename,
                     std::string description = "grid file") const;

  /// Write the grid data without labels, as they are represented in memory
  std::ostream & write_opendx(std::ostream &os) const;

  /// Write the grid data without labels, as they are represented in memory
  int write_opendx(std::string const &filename,
                   std::string description = "grid file") const;
};


enum class scalar_type {
  real,
  integer
};

template <scalar_type T>
struct scalar_type_map;

template <>
struct scalar_type_map<scalar_type::real> {
  using type = cvm::real;
};

template <>
struct scalar_type_map<scalar_type::integer> {
  using type = size_t;
};

template <scalar_type T>
class colvar_grid_scalar_function : public colvar_grid<typename scalar_type_map<T>::type>
{
public:
  using base_class = colvar_grid<typename scalar_type_map<T>::type>;
  using scalar_type_class = typename scalar_type_map<T>::type;
  /// Default constructor
  colvar_grid_scalar_function(): base_class(){};

  /// Copy constructor (needed because of the grad pointer)
  colvar_grid_scalar_function(colvar_grid_scalar_function const &g) : base_class(g){};
  colvar_grid_scalar_function(std::vector<colvar *>  &colvars,std::string config): base_class(colvars, 0, 1, false, nullptr,
    config){};

  /// Destructor
  virtual ~colvar_grid_scalar_function(){};

  /// Constructor from a vector of colvars
  colvar_grid_scalar_function(std::vector<colvar *> &colvars,
                     std::shared_ptr<const colvar_grid_params> params = nullptr,
                     bool add_extra_bin = false,
                     std::string config = std::string()): base_class(colvars, 0.0, 1, add_extra_bin, params, config){};

  /// Constructor from a multicol file
  colvar_grid_scalar_function(std::string const &filename):base_class(filename, 1){};

  /// Accumulate the value
  inline void acc_value(std::vector<int> const &ix,
                        cvm::real const &new_value,
                        size_t const &imult = 0)
  {
    (void) imult;
    // only legal value of imult here is 0
    this->data[this->address(ix)] += new_value;
    this->has_data = true;
  }

  /// \brief Return the gradient of the scalar field from finite differences
  /// Input coordinates are those of gradient grid, shifted wrt scalar grid
  /// Should not be called on edges of scalar grid, provided the latter has
  /// margins (extra bins) wrt gradient grid
  inline void vector_gradient_finite_diff( const std::vector<int> &ix0, std::vector<cvm::real> &grad)
  {
    cvm::real A0, A1;
    std::vector<int> ix;
    size_t i, j, k, n;

     if (this->nd == 2) {
      for (n = 0; n < 2; n++) {
        ix = ix0;
        A0 = this->value(ix);
        ix[n]++; this->wrap(ix);
        A1 = this->value(ix);
        ix[1-n]++; this->wrap(ix);
        A1 += this->value(ix);
        ix[n]--; this->wrap(ix);
        A0 += this->value(ix);
        grad[n] = 0.5 * (A1 - A0) / this->widths[n];
      }
    } else if (this->nd == 3) {

      cvm::real p[8]; // potential values within cube, indexed in binary (4 i + 2 j + k)
      ix = ix0;
      int index = 0;
      for (i = 0; i<2; i++) {
        ix[1] = ix0[1];
        for (j = 0; j<2; j++) {
          ix[2] = ix0[2];
          for (k = 0; k<2; k++) {
            this->wrap(ix);
            p[index++] = this->value(ix);
            ix[2]++;
          }
          ix[1]++;
        }
        ix[0]++;
      }

      // The following would be easier to read using binary literals
      //                  100    101    110    111      000    001    010   011
      grad[0] = 0.25 * ((p[4] + p[5] + p[6] + p[7]) - (p[0] + p[1] + p[2] + p[3])) / this->widths[0];
      //                  010     011    110   111      000    001    100   101
      grad[1] = 0.25 * ((p[2] + p[3] + p[6] + p[7]) - (p[0] + p[1] + p[4] + p[5])) / this->widths[1];
      //                  001    011     101   111      000    010   100    110
      grad[2] = 0.25 * ((p[1] + p[3] + p[5] + p[7]) - (p[0] + p[2] + p[4] + p[6])) / this->widths[2];
    } else {
      cvm::error("Finite differences available in dimension 2 and 3 only.");
    }
  }


  /// \brief Return the log-gradient from finite differences
  /// on the *same* grid for dimension n
  /// (colvar_grid_scalar)
  inline cvm::real log_gradient_finite_diff(const std::vector<int> &ix0,
                                            int n = 0, int offset = 0)
  {
    cvm::real A0, A1, A2;
    std::vector<int> ix = ix0;

    // TODO this can be rewritten more concisely with wrap_edge()
    if (this->periodic[n]) {
      ix[n]--; this->wrap(ix);
      A0 = this->value(ix) + offset;
      ix = ix0;
      ix[n]++; this->wrap(ix);
      A1 = this->value(ix) + offset;
      if (A0 * A1 == 0) {
        return 0.; // can't handle empty bins
      } else {
        return (cvm::logn(A1) - cvm::logn(A0))
          / (this->widths[n] * 2.);
      }
    } else if (ix[n] > 0 && ix[n] < this->nx[n]-1) { // not an edge
      ix[n]--;
      A0 = this->value(ix) + offset;
      ix = ix0;
      ix[n]++;
      A1 = this->value(ix) + offset;
      if (A0 * A1 == 0) {
        return 0.; // can't handle empty bins
      } else {
        return (cvm::logn(A1) - cvm::logn(A0))
          / (this->widths[n] * 2.);
      }
    } else {
      // edge: use 2nd order derivative
      int increment = (ix[n] == 0 ? 1 : -1);
      // move right from left edge, or the other way around
      A0 = this->value(ix) + offset;
      ix[n] += increment; A1 = this->value(ix) + offset;
      ix[n] += increment; A2 = this->value(ix) + offset;
      if (A0 * A1 * A2 == 0) {
        return 0.; // can't handle empty bins
      } else {
        return (-1.5 * cvm::logn(A0) + 2. * cvm::logn(A1)
          - 0.5 * cvm::logn(A2)) * increment / this->widths[n];
      }
    }
  }


  /// \brief Return the gradient of discrete count from finite differences
  /// on the *same* grid for dimension n
  /// (colvar_grid_scalar)
  inline cvm::real gradient_finite_diff(const std::vector<int> &ix0,
                                        int n = 0)
  {
    cvm::real Am1, A0, A1, A2;
    std::vector<int> ix = ix0;

    // FIXME this can be rewritten more concisely with wrap_edge()
    if (this->periodic[n]) {
      ix[n]--; this->wrap(ix);
      Am1 = this->value(ix);
      ix = ix0; this->wrap(ix);
      A0 = this->value(ix); // Used just to detect gaps
      ix[n]++; this->wrap(ix);
      A1 = this->value(ix);
      if (Am1 * A0 * A1 == 0) {
        return 0.; // can't handle empty bins
      } else {
        return (A1 - Am1) / (this->widths[n] * 2.);
      }
    } else if (ix[n] > 0 && ix[n] < this->nx[n]-1) { // not an edge
      ix[n]--;
      Am1 = this->value(ix);
      ix = ix0;
      A0 = this->value(ix);
      ix[n]++;
      A1 = this->value(ix);
      if (Am1 * A0 * A1 == 0) {
        return 0.; // can't handle empty bins
      } else {
        return (A1 - Am1) / (this->widths[n] * 2.);
      }
    } else {
      // edge: use 2nd order derivative
      int increment = (ix[n] == 0 ? 1 : -1);
      // move right from left edge, or the other way around
      A0 = this->value(ix);
      ix[n] += increment; A1 = this->value(ix);
      ix[n] += increment; A2 = this->value(ix);
      return (-1.5 * A0 + 2. * A1
          - 0.5 * A2) * increment / this->widths[n];
    }
  }

  using base_class::value_output;
  /// \brief Return the value of the function at ix divided by its
  /// number of samples (if the count grid is defined)
  inline scalar_type_class value_output(std::vector<int> const &ix,
                                        size_t const &imult = 0) const
  {
    if (imult > 0) {
      cvm::error("Error: trying to access a component "
                 "larger than 1 in a scalar data grid.\n");
      return 0.;
    }
    else {
      return this->data[this->address(ix) + imult];
    }
  }
  using base_class::value_input;
  /// Enter or add value but also deal with count grid
  void value_input(std::vector<int> const &ix,
                           cvm::real const &new_value,
                           size_t const &imult = 0,
                           bool add = false)
  {
    if (imult > 0) {
      cvm::error("Error: trying to access a component "
                 "larger than 1 in a scalar data grid.\n");
      return;
    }
    if (add) {
        this->data[this->address(ix)] += new_value;
        if (this->has_parent_data) {
          // save newly read data for inputting parent grid
          this->new_data[this->address(ix)] = static_cast<scalar_type_class>(new_value);
        }
    } else {
        this->data[this->address(ix)] = new_value;
    }
    this->has_data = true;
  }


  /// \brief Return the highest value
  cvm::real maximum_value() const{
    cvm::real max = this->data[0];
    for (size_t i = 0; i < this->nt; i++) {
      if (this->data[i] > max)
        max = this->data[i];
    }
    return max;
  };

  /// \brief Return the lowest value
  cvm::real minimum_value() const{
    cvm::real min = this->data[0];
    for (size_t i = 0; i < this->nt; i++) {
      if (this->data[i] < min)
        min = this->data[i];
    }
    return min;
  };

  /// \brief Return the lowest positive value
  cvm::real minimum_pos_value() const{
    cvm::real minpos = this->data[0];
    size_t i;
    for (i = 0; i < this->nt; i++) {
      if (this->data[i] > 0) {
        minpos = this->data[i];
        break;
      }
    }
    for (i = 0; i < this->nt; i++) {
      if (this->data[i] > 0 && this->data[i] < minpos)
        minpos = this->data[i];
    }
    return minpos;
  };

  /// \brief Calculates the integral of the map (uses widths if they are defined)
  cvm::real integral() const{
    cvm::real sum = 0.0;
    for (size_t i = 0; i < this->nt; i++) {
      sum += this->data[i];
    }
    cvm::real bin_volume = 1.0;
    for (size_t id = 0; id < this->widths.size(); id++) {
      bin_volume *= this->widths[id];
    }
    return bin_volume * sum;
  };

  /// \brief Assuming that the map is a normalized probability density,
  ///        calculates the entropy (uses widths if they are defined)
  cvm::real entropy() const{
    cvm::real sum = 0.0;
    for (size_t i = 0; i < this->nt; i++) {
      if (this->data[i] > 0) {
        sum += -1.0 * this->data[i] * cvm::logn(this->data[i]);
      }
    }
    cvm::real bin_volume = 1.0;
    for (size_t id = 0; id < this->widths.size(); id++) {
      bin_volume *= this->widths[id];
    }
    return bin_volume * sum;
  };

  /// \brief Return the average number of samples in a given "radius" around current bin
  /// Really a hypercube of length 2*radius + 1
  inline int local_sample_count(int radius)
  {
    std::vector<int> ix0 = this->new_index();
    std::vector<int> ix = this->new_index();

    for (size_t i = 0; i < this->nd; i++) {
      ix0[i] = this->current_bin_scalar_bound(i);
    }
    if (radius < 1) {
      // Simple case: no averaging
      if (this->index_ok(ix0))
        return this->value(ix0);
      else
        return 0;
    }
    scalar_type_class count = 0;
    size_t nbins = 0;
    int i, j, k;
    bool edge;
    ix = ix0;
    // Treat each dimension separately to simplify code
    switch (this->nd)
    {
      case 1:
        for (i = -radius; i <= radius; i++) {
          ix[0] = ix0[0] + i;
          edge = this->wrap_detect_edge(ix);
          if (!edge) {
            nbins++;
            count += this->value(ix);
          }
        }
        break;
      case 2:
        for (i = -radius; i <= radius; i++) {
          ix[0] = ix0[0] + i;
          for (j = -radius; j <= radius; j++) {
            ix[1] = ix0[1] + j;
            edge = this->wrap_detect_edge(ix);
            if (!edge) {
              nbins++;
              count += this->value(ix);
            }
          }
        }
        break;
      case 3:
        for (i = -radius; i <= radius; i++) {
          ix[0] = ix0[0] + i;
          for (j = -radius; j <= radius; j++) {
            ix[1] = ix0[1] + j;
            for (k = -radius; k <= radius; k++) {
              ix[2] = ix0[2] + k;
              edge = this->wrap_detect_edge(ix);
              if (!edge) {
                nbins++;
                count += this->value(ix);
              }
            }
          }
        }
        break;
      default:
        cvm::error("Error: local_sample_count is not implemented for grids of dimension > 3", COLVARS_NOT_IMPLEMENTED);
        break;
    }

    if (nbins)
      // Integer division - an error on the order of 1 doesn't matter
        return count / nbins;
    else
      return 0.0;
  }


  /// \brief Return the RMSD between this grid's data and another one
  /// Grids must have the same dimensions.
  cvm::real grid_rmsd(colvar_grid_scalar_function const &other_grid) const{
    if (other_grid.data.size() != this->data.size()) {
      cvm::error("Error: trying to subtract two grids with "
                 "different size.\n");
      return -1.;
    }

    cvm::real sum2 = 0.0;
    for (size_t i = 0; i < this->data.size(); i++) {
      cvm::real d = other_grid.data[i] - this->data[i];
      sum2 += d * d;
    }
    return sqrt(sum2 / this->data.size());
  };


  void increase(std::vector<int> const &ix, scalar_type_class fact = 1)
  {
    this->data[this->address(ix)] += fact;
  }
  /// \brief Get the binned count indexed by ix from the newly read data
  inline scalar_type_class const & new_value(std::vector<int> const &ix)
  {
    return this->new_data[this->address(ix)];
  }
};


typedef colvar_grid_scalar_function<scalar_type::real> colvar_grid_scalar;
typedef colvar_grid_scalar_function<scalar_type::integer> colvar_grid_count;

/// Class for accumulating the gradient of a scalar function on a grid
class colvar_grid_gradient : public colvar_grid<cvm::real>
{
public:

  /// \brief Provide the sample count by which each binned value
  /// should be divided
  std::shared_ptr<colvar_grid_count> samples;
  /// Continuous weights instead of samples, eg. for smoothed ABF
  std::shared_ptr<colvar_grid_scalar> weights;
  /// Default constructor
  colvar_grid_gradient();

  /// Destructor
  virtual ~colvar_grid_gradient()
  {}

  // /// Constructor from specific sizes arrays
  // colvar_grid_gradient(std::vector<int> const &nx_i);

  // /// Constructor from a vector of colvars
  // colvar_grid_gradient(std::vector<colvar *>  &colvars,
  //                      std::string config = std::string());

  /// Constructor from a multicol file
  explicit colvar_grid_gradient(std::string const &filename, std::shared_ptr<colvar_grid_count> samples_in = nullptr);
  explicit colvar_grid_gradient(std::string const &filename, std::shared_ptr<colvar_grid_scalar> weights_in= nullptr);


  /// Constructor from a vector of colvars and a pointer to the count grid
  explicit colvar_grid_gradient(std::vector<colvar *> &colvars,
                       std::shared_ptr<colvar_grid_count> samples_in = nullptr,
                       std::shared_ptr<const colvar_grid_params> params = nullptr,
                       std::string config = std::string());
  explicit  colvar_grid_gradient(std::vector<colvar *> &colvars,
                       std::shared_ptr<colvar_grid_scalar> weights_in,
                       std::shared_ptr<const colvar_grid_params> params = nullptr,
                       std::string config = std::string());


  /// Parameters for smoothing data with low sampling
  int full_samples;
  int min_samples;
  // TOD: maybe make those class members
  // int smoothing;
  // int cutoff;

  /// Write the current grid parameters to a string
  std::string get_state_params() const override;

  /// Read new grid parameters from a string
  int parse_params(std::string const &conf,
                   colvarparse::Parse_Mode const parse_mode = colvarparse::parse_normal);

  /// Read all grid parameters and data from a formatted stream
  std::istream & read_restart(std::istream &is);

  /// Read all grid parameters and data from an unformatted stream
  cvm::memory_stream & read_restart(cvm::memory_stream &is);

  /// Write all grid parameters and data to a formatted stream
  std::ostream & write_restart(std::ostream &os);

  /// Write all grid parameters and data to an unformatted stream
  cvm::memory_stream & write_restart(cvm::memory_stream &os);

  /// Read all grid parameters and data from a formatted stream
  std::istream &read_raw(std::istream &is);

  /// Read all grid parameters and data from an unformatted stream
  cvm::memory_stream &read_raw(cvm::memory_stream &is);

  /// Write all grid data to a formatted stream (without labels, as they are represented in memory)
  /// \param[in,out] os Stream object
  /// \param[in] buf_size Number of values per line
  std::ostream &write_raw(std::ostream &os, size_t const buf_size = 3) const;

  /// Write all grid data to an unformatted stream
  /// \param[in,out] os Stream object
  /// \param[in] buf_size Number of values per line (note: ignored because there is no formatting)
  cvm::memory_stream &write_raw(cvm::memory_stream &os, size_t const buf_size = 3) const;

  /// Read a grid written by write_multicol(), incrementin if data is true
  std::istream & read_multicol(std::istream &is, bool add = false);

  /// Read a grid written by write_multicol(), incrementin if data is true
  int read_multicol(std::string const &filename,
                            std::string description = "grid file",
                            bool add = false);

  /// Write grid in a format which is both human-readable and gnuplot-friendly
  std::ostream & write_multicol(std::ostream &os) const;

  /// Write grid in a format which is both human-readable and gnuplot-friendly
  int write_multicol(std::string const &filename,
                             std::string description = "grid file") const;

  /// Write the grid data without labels, as they are represented in memory
  std::ostream & write_opendx(std::ostream &os) const;

  /// Write the grid data without labels, as they are represented in memory
  int write_opendx(std::string const &filename,
                           std::string description = "grid file") const;

  /// \brief Get a vector with the binned value(s) indexed by ix, normalized if applicable
  inline void vector_value(std::vector<int> const &ix, std::vector<cvm::real> &v) const
  {
    cvm::real const * p = &value(ix);

    if (samples) {
      int count = samples->value(ix);
      if (count) {
        cvm::real invcount = 1.0 / count;
        for (size_t i = 0; i < mult; i++) {
          v[i] = invcount * p[i];
        }
      } else {
        for (size_t i = 0; i < mult; i++) {
          v[i] = 0.0;
        }
      }
    }
    else if (weights) {
      cvm::real weight = weights->value(ix);
      if (weight > 0) {
        // cvm::log( "we use weights" + cvm::to_str(weight) );
        // cvm::log("gradients is : [" + cvm::to_str(p[0]) + " " + cvm::to_str(p[1]) + "]");
        cvm::real invcount = 1.0 / weight;
        for (size_t i = 0; i < mult; i++) {
          v[i] = invcount * p[i];
        }
        // cvm::log("gradients is : [" + cvm::to_str(v[0]) + " " + cvm::to_str(v[1]) + "]");
      } else {
        for (size_t i = 0; i < mult; i++) {
          v[i] = 0.0;
        }
      }
    } else {
      for (size_t i = 0; i < mult; i++) {
        v[i] = p[i];
      }
    }
  }


  /// \brief Accumulate the value
  inline void acc_value(std::vector<int> const &ix, std::vector<colvarvalue> const &values) {
    for (size_t imult = 0; imult < mult; imult++) {
      data[address(ix) + imult] += values[imult].real_value;
    }
    if (samples) {
      samples->increase(ix,1);
    }
    else if (weights) {
      cvm::error("acc_value is only in use for ti bias and shouldn't be used with non-integer weights", COLVARS_ERROR);
    }
  }

  /// \brief Accumulate the gradient based on the force (i.e. sums the
  /// opposite of the force)
  inline void acc_force(std::vector<int> const &ix, cvm::real const *forces, bool b_smoothed = false,
                              cvm::real fact = 1.0) {
    for (size_t imult = 0; imult < mult; imult++) {
      data[address(ix) + imult] -= forces[imult] * fact;
    }
    if (samples) {
      samples->increase(ix, 1);
    }
    else if (weights) {
      weights->increase(ix,fact);
    }
    // cvm::log("we increase count to " + cvm::to_str(samples ? samples->value(ix) : weights->value(ix)));
    // cvm::log("force input : " + cvm::to_str(*forces) + " reweighted forces = " + cvm::to_str(std::vector<cvm::real>{data[address(ix)], data[address(ix) + 1]}));
  }
  bool index_ok_two_vec(std::vector<int> const &ix, std::vector<int> const ix_min, std::vector<int> const ix_max) const
  {
    for (size_t i = 0; i < nd; i++) {
      if ( (ix[i] < ix_min[i]) || (ix[i] > ix_max[i]))
        return false;
    }
    return true;
  }
  inline void incr_two_vec(std::vector<int> &ix, std::vector<int> const ix_min, std::vector<int> const ix_max) const
  {
    for (int i = ix.size()-1; i >= 0; i--) {

      ix[i]++;

      if (ix[i] > ix_max[i]) {

        if (i > 0) {
          ix[i] = ix_min[i];
          continue;
        } else {
          // this is the last iteration, a non-valid index is being
          // set for the outer index, which will be caught by
          // index_ok()
          // ix[0] = ix_min[0];
          return;
        }
      } else {
        return;
      }
    }
  }
    /// \brief Accumulate the gradient based on the force (i.e. sums the
  /// opposite of the force)
  /// TODO: explain what CV value is
  inline void acc_abf_force(std::vector<cvm::real> const &cv_value, std::vector<int> const &bin_value,
                        cvm::real const *force,
                        bool b_smoothed,
                        cvm::real smoothing = 2.0) {



    if (b_smoothed && weights->value(bin_value) < full_samples) {
      cvm::real ixmin, ixmax;
      std::vector<int> bin(nd, 0);
      // We process points where exp(-d^2 / 2sigma^2) > 10^-3
      // This implies distance < 3.72 * sigma
      cvm::real inv_squared_smooth = 1/ (std::max(smoothing*smoothing, 1e-5));
      cvm::real const cutoff_factor = 3.72;
      // TODO: make sure that this is not > min nx /2
      cvm::real const cutoff = cutoff_factor * smoothing; // take like floor()
      std::vector<int>ix_min(nd, 0);
      std::vector<int>ix_max(nd, 0);
      std::vector<int>periodic_offset(nd,0);
      for (size_t i =0; i < nd; i++) {
        ixmin = cv_value[i] - cutoff;
        ixmax = cv_value[i] + cutoff;
        int cut_ixmin = static_cast<int>(cvm::floor(std::max(ixmin,static_cast<cvm::real>(0))));
        int cut_ixmax = static_cast<int>(cvm::floor(std::min(ixmax, static_cast<cvm::real>(nx[0]) - 1)));
        if (periodic[i]) {
          periodic_offset[i] = cut_ixmin - static_cast<int>(cvm::floor(ixmin)) - static_cast<int>(cvm::floor(ixmax)) + cut_ixmax;
          cut_ixmin = cut_ixmin + periodic_offset[i];
          cut_ixmax = cut_ixmax + periodic_offset[i];
        }
        ix_min[i] = cut_ixmin;
        ix_max[i] = cut_ixmax;
      }
      cvm::real dist = 0;
      std::vector<bool>stop_vec(nd, false);
      // TODO: we can parallelize here i think
      for (std::vector<int> ix = ix_min; index_ok_two_vec(ix, ix_min, ix_max); incr_two_vec(ix, ix_min, ix_max)) {
        dist = 0;
        std::vector<int>ix_copy = std::vector<int>(ix);
        for (size_t dim = 0; dim < nd; dim++) {
            if (periodic[dim]) {
              ix_copy[dim] = (ix_copy[dim] - periodic_offset[dim]) % nx[dim];
            }
          dist += std::pow(ix_copy[dim] +0.5 - cv_value[dim], 2);
        }
        acc_force(ix_copy, force, true, cvm::exp(-dist * inv_squared_smooth));
      }
    } else {
      // cvm::log("we update the force with acc_force, force = " + cvm::to_str(force[0]) + " " + cvm::to_str(force[1]));
      acc_force(bin_value, force);
    }
  }

  /// \brief Return the value of the function at ix divided by its
  /// number of samples (if the count grid is defined)
  virtual cvm::real value_output(std::vector<int> const &ix,
                                 size_t const &imult = 0) const override
  {
    cvm::real s;
    if (samples || weights) {
      s = samples ? static_cast<cvm::real>(samples->value(ix)) : weights->value(ix);
      return s > 0 ?(data[address(ix) + imult] / s) : 0.0;
    }
    return data[address(ix) + imult];
  }

  /// Compute the inverse weight corresponding to smoothing factor as in ABF
  /// to normalize sums over steps into averages
  inline cvm::real smooth_inverse_weight(cvm::real weight)
  {
    cvm::real fact;
    if ( weight <= min_samples ) {
      fact = 0.0;
    } else if ( weight < full_samples ) {
      fact = (weight - min_samples) / (weight * cvm::real(full_samples - min_samples));
    } else {
      fact = 1.0 / weight;
    }
    return fact;
  }


  /// \brief Return the scalar value of the function at ix divided by its
  /// number of samples (if the count grid is defined), possibly smoothed
  /// by a ramp function going from 0 to 1 between minSamples and fullSamples.
  /// Only makes sense if dimension is 1
  virtual inline cvm::real value_output_smoothed(std::vector<int> const &ix, bool smoothed = true)
  {
    cvm::real weight, fact;

    if (samples || weights) {
      weight = samples ? cvm::real(samples->value(ix)) : weights->value(ix);
    } else {
      weight = 1.;
    }

    if (smoothed) {
      fact = smooth_inverse_weight(weight);
    } else {
      fact = weight > 0. ? 1. / weight : 0.;
    }

    return fact * data[address(ix)];
  }


  /// \brief Obtain the vector value of the function at ix divided by its
  /// number of samples (if the count grid is defined), possibly smoothed
  /// by a ramp function going from 0 to 1 between minSamples and fullSamples.
  inline void vector_value_smoothed(std::vector<int> const &ix, cvm::real *grad, bool smoothed = true)
  {
    cvm::real weight, fact;

    if (samples) {
      weight = static_cast<cvm::real>(samples->value(ix));
    } else if (weights) {
      weight = weights->value(ix);
    } else {
      weight = 1.;
    }
    if (smoothed) {
      fact = smooth_inverse_weight(weight);
    } else {
      fact = weight > 0. ? 1. / weight : 0.;
    }
    cvm::real *p = &(data[address(ix)]);

    for (size_t imult = 0; imult < mult; imult++) {
      grad[imult] = fact * p[imult];
    }
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (it may have been rescaled or
  /// manipulated)
  virtual void value_input(std::vector<int> const &ix,
                           cvm::real const &new_value,
                           size_t const &imult = 0,
                           bool add = false) override
  {
    if (add) {
      if (samples || weights) {
        data[address(ix) + imult] += new_value * (samples ? samples->new_value(ix) : weights->new_value(ix));
      }
      else
        data[address(ix) + imult] += new_value;
    } else {
      if (samples || weights)
        data[address(ix) + imult] = new_value * (samples ? samples->value(ix) : weights->value(ix));
      else
        data[address(ix) + imult] = new_value;
    }
    has_data = true;
  }


  /// Compute and return average value for a 1D gradient grid
  inline cvm::real average(bool smoothed = false)
  {
    if (nd != 1 || nx[0] == 0) {
      return 0.0;
    }

    cvm::real sum = 0.0;
    for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix)) {
      sum += value_output_smoothed(ix, smoothed);
    }

    return (sum / cvm::real(nx[0]));
  }

  /// \brief Return the RMSD between this grid's data and another one
  /// Grids must have the same dimensions.
  cvm::real grid_rmsd(colvar_grid_gradient const &other_grid) const;

  /// \brief If the grid is 1-dimensional, integrate it and write the
  /// integral to a file (DEPRECATED by the colvargrid_integrate class)
  void write_1D_integral(std::ostream &os);

};
 #include "colvargrid_def.h"

 #endif
