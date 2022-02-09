// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARBIAS_H
#define COLVARBIAS_H

#include "colvar.h"
#include "colvarparse.h"
#include "colvardeps.h"

class colvar_grid_scalar;

/// \brief Collective variable bias, base class
class colvarbias
  : public virtual colvarparse, public virtual colvardeps {
public:

  /// Name of this bias
  std::string name;

  /// Keyword indicating the type of this bias
  std::string bias_type;

  /// Keyword used in state files (== bias_type most of the time)
  std::string state_keyword;

  /// Track how many times a bias of this type was defined
  int rank;

  /// Add a new collective variable to this bias
  int add_colvar(std::string const &cv_name);

  /// How many variables are defined for this bias
  inline size_t num_variables() const
  {
    return colvars.size();
  }

  /// Access the variables vector
  inline std::vector<colvar *> *variables()
  {
    return &colvars;
  }

  /// Access the i-th variable
  inline colvar * variables(int i) const
  {
    return colvars[i];
  }

  /// Retrieve colvar values and calculate their biasing forces
  /// Some implementations may use calc_energy() and calc_forces()
  virtual int update();

  /// Returns true if the current step represent a valid increment, whose data
  /// can be recorded (as opposed to e.g. a continuation step from a restart)
  virtual bool can_accumulate_data();

  /// Compute the energy of the bias
  /// Uses the vector of colvar values provided if not NULL, and the values
  /// currently cached in the bias instance otherwise
  virtual int calc_energy(std::vector<colvarvalue> const *values);

  /// Compute the forces due to the bias
  /// Uses the vector of colvar values provided if not NULL, and the values
  /// currently cached in the bias instance otherwise
  virtual int calc_forces(std::vector<colvarvalue> const *values);

  /// Send forces to the collective variables
  void communicate_forces();

  /// Carry out operations needed before next step is run
  virtual int end_of_step();

  /// Load new configuration - force constant and/or centers only
  virtual int change_configuration(std::string const &conf);

  /// Calculate change in energy from using alternate configuration
  virtual cvm::real energy_difference(std::string const &conf);

  /// Give the total number of bins for a given bias.
  // FIXME this is currently 1D only
  virtual int bin_num();
  /// Calculate the bin index for a given bias.
  // FIXME this is currently 1D only
  virtual int current_bin();
  //// Give the count at a given bin index.
  // FIXME this is currently 1D only
  virtual int bin_count(int bin_index);
  //// Share information between replicas, whatever it may be.
  virtual int replica_share();

  /// Perform analysis tasks
  virtual void analyze() {}

  /// \brief Constructor
  colvarbias(char const *key);

  /// \brief Parse config string and (re)initialize
  virtual int init(std::string const &conf);

  /// \brief Initialize dependency tree
  virtual int init_dependencies();

  /// \brief Set to zero all mutable data
  virtual int reset();

private:

  /// Default constructor
  colvarbias();

  /// Copy constructor
  colvarbias(colvarbias &);

public:

  /// \brief Delete everything
  virtual int clear();

  /// \brief Delete only the allocatable data (save memory)
  virtual int clear_state_data();

  /// Destructor
  virtual ~colvarbias();

  /// Write the values of specific mutable properties to a string
  virtual std::string const get_state_params() const;

  /// Read the values of specific mutable properties from a string
  virtual int set_state_params(std::string const &state_conf);

  /// Write all mutable data not already written by get_state_params()
  virtual std::ostream & write_state_data(std::ostream &os)
  {
    return os;
  }

  /// Read all mutable data not already set by set_state_params()
  virtual std::istream & read_state_data(std::istream &is)
  {
    return is;
  }

  /// Read a keyword from the state data (typically a header)
  /// \param Input stream
  /// \param Keyword labeling the header block
  std::istream & read_state_data_key(std::istream &is, char const *key);

  /// Write the bias configuration to a state file or other stream
  std::ostream & write_state(std::ostream &os);

  /// Read the bias configuration from a restart file or other stream
  std::istream & read_state(std::istream &is);

  /// Write the bias state to a file with the given prefix
  int write_state_prefix(std::string const &prefix);

  /// Write the bias state to a string
  int write_state_string(std::string &output);

  /// Read the bias state from a file with this name or prefix
  int read_state_prefix(std::string const &prefix);

  /// Read the bias state from this string buffer
  int read_state_string(char const *buffer);

  /// Write a label to the trajectory file (comment line)
  virtual std::ostream & write_traj_label(std::ostream &os);

  /// Output quantities such as the bias energy to the trajectory file
  virtual std::ostream & write_traj(std::ostream &os);

  /// (Re)initialize the output files (does not write them yet)
  virtual int setup_output()
  {
    return COLVARS_OK;
  }

  /// Frequency for writing output files
  size_t output_freq;

  /// Write any output files that this bias may have (e.g. PMF files)
  virtual int write_output_files()
  {
    return COLVARS_OK;
  }

  /// Use this prefix for all output files
  std::string output_prefix;

  /// If this bias is communicating with other replicas through files, send it to them
  virtual int write_state_to_replicas()
  {
    return COLVARS_OK;
  }

  inline cvm::real get_energy()
  {
    return bias_energy;
  }

  /// \brief Implementation of the feature list for colvarbias
  static std::vector<feature *> cvb_features;

  /// \brief Implementation of the feature list accessor for colvarbias
  virtual const std::vector<feature *> &features() const
  {
    return cvb_features;
  }
  virtual std::vector<feature *> &modify_features()
  {
    return cvb_features;
  }
  static void delete_features() {
    for (size_t i=0; i < cvb_features.size(); i++) {
      delete cvb_features[i];
    }
    cvb_features.clear();
  }

protected:

  /// \brief Pointers to collective variables to which the bias is
  /// applied; current values and metric functions will be obtained
  /// through each colvar object
  std::vector<colvar *>    colvars;

  /// \brief Up to date value of each colvar
  std::vector<colvarvalue> colvar_values;

  /// \brief Current forces from this bias to the variables
  std::vector<colvarvalue> colvar_forces;

  /// \brief Forces last applied by this bias to the variables
  std::vector<colvarvalue> previous_colvar_forces;

  /// \brief Current energy of this bias (colvar_forces should be obtained by deriving this)
  cvm::real                bias_energy;

  /// Whether to write the current bias energy from this bias to the trajectory file
  bool                     b_output_energy;

  /// \brief Whether this bias has already accumulated information
  /// (for history-dependent biases)
  bool                     has_data;

  /// \brief Step number read from the last state file
  cvm::step_number         state_file_step;

  /// Flag used to tell if the state string being read is for this bias
  bool matching_state;

  /// \brief The biasing forces will be scaled by the factor in this grid
  /// if b_bias_force_scaled is true
  colvar_grid_scalar*      biasing_force_scaling_factors;
  std::vector<int>         biasing_force_scaling_factors_bin;
};


class colvar_grid_gradient;
class colvar_grid_count;

/// \brief Base class for unconstrained thermodynamic-integration FE estimator
class colvarbias_ti : public virtual colvarbias {
public:

  colvarbias_ti(char const *key);
  virtual ~colvarbias_ti();

  virtual int clear_state_data();

  virtual int init(std::string const &conf);
  virtual int init_grids();
  virtual int update();

  /// Subtract applied forces (either last forces or argument) from the total
  /// forces
  virtual int update_system_forces(std::vector<colvarvalue> const
                                   *subtract_forces);

  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &state_conf);
  virtual std::ostream & write_state_data(std::ostream &os);
  virtual std::istream & read_state_data(std::istream &is);
  virtual int write_output_files();

protected:

  /// \brief Forces exerted from the system to the associated variables
  std::vector<colvarvalue> ti_system_forces;

  /// Averaged system forces
  colvar_grid_gradient *ti_avg_forces;

  /// Histogram of sampled data
  colvar_grid_count *ti_count;

  /// Because total forces may be from the last simulation step,
  /// store the index of the variables then
  std::vector<int> ti_bin;
};

#endif
