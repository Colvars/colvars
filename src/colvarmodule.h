// -*- c++ -*-

#ifndef COLVARMODULE_H
#define COLVARMODULE_H

#include "colvars_version.h"

#ifndef COLVARS_DEBUG
#define COLVARS_DEBUG false
#endif

/*! \mainpage Main page
This is the Developer's documentation for the Collective Variables Module.

You can browse the class hierarchy or the list of source files.
 */

/// \file colvarmodule.h
/// \brief Collective variables main module
///
/// This file declares the main class for defining and manipulating
/// collective variables: there should be only one instance of this
/// class, because several variables are made static (i.e. they are
/// shared between all object instances) to be accessed from other
/// objects.

#define COLVARS_OK 0
#define COLVARS_ERROR   1
#define COLVARS_NOT_IMPLEMENTED (1<<1)
#define INPUT_ERROR     (1<<2) // out of bounds or inconsistent input
#define BUG_ERROR       (1<<3) // Inconsistent state indicating bug
#define FILE_ERROR      (1<<4)
#define MEMORY_ERROR    (1<<5)
#define FATAL_ERROR     (1<<6) // Should be set, or not, together with other bits
#define DELETE_COLVARS  (1<<7) // Instruct the caller to delete cvm
#define COLVARS_NO_SUCH_FRAME (1<<8) // Cannot load the requested frame

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <list>

#ifdef NAMD_VERSION
// use Lustre-friendly wrapper to POSIX write()
#include "fstream_namd.h"
#endif

class colvarparse;
class colvar;
class colvarbias;
class colvarproxy;
class colvarscript;


/// \brief Collective variables module (main class)
///
/// Class to control the collective variables calculation.  An object
/// (usually one) of this class is spawned from the MD program,
/// containing all i/o routines and general interface.
///
/// At initialization, the colvarmodule object creates a proxy object
/// to provide a transparent interface between the MD program and the
/// child objects
class colvarmodule {

private:

  /// Impossible to initialize the main object without arguments
  colvarmodule();

public:

  friend class colvarproxy;
  // TODO colvarscript should be unaware of colvarmodule's internals
  friend class colvarscript;

  /// Defining an abstract real number allows to switch precision
  typedef  double    real;
  /// Residue identifier
  typedef  int       residue_id;

  class rvector;
  template <class T> class vector1d;
  template <class T> class matrix2d;
  class quaternion;
  class rotation;

  /// \brief Atom position (different type name from rvector, to make
  /// possible future PBC-transparent implementations)
  typedef rvector atom_pos;

  /// \brief 3x3 matrix of real numbers
  class rmatrix;

  // allow these classes to access protected data
  class atom;
  class atom_group;
  friend class atom;
  friend class atom_group;
  typedef std::vector<atom>::iterator       atom_iter;
  typedef std::vector<atom>::const_iterator atom_const_iter;

  /// Module-wide error state
  /// see constants at the top of this file
protected:

  static int errorCode;

public:

  static void set_error_bits(int code);

  static bool get_error_bit(int code);

  static inline int get_error()
  {
    return errorCode;
  }

  static void clear_error();


  /// Current step number
  static long it;
  /// Starting step number for this run
  static long it_restart;

  /// Return the current step number from the beginning of this run
  static inline long step_relative()
  {
    return it - it_restart;
  }

  /// Return the current step number from the beginning of the whole
  /// calculation
  static inline long step_absolute()
  {
    return it;
  }

  /// If true, get it_restart from the state file; if set to false,
  /// the MD program is providing it
  bool it_restart_from_state_file;

  /// \brief Finite difference step size (if there is no dynamics, or
  /// if gradients need to be tested independently from the size of
  /// dt)
  static real debug_gradients_step_size;

private:

  /// Prefix for all output files for this run
  std::string cvm_output_prefix;

public:
  /// Accessor for the above
  static inline std::string &output_prefix()
  {
    colvarmodule *cv = colvarmodule::main();
    return cv->cvm_output_prefix;
  }

private:

  /// Array of collective variables
  std::vector<colvar *> colvars;

  /// Array of collective variables
  std::vector<colvar *> colvars_active;

  /// Collective variables to be calculated on different threads;
  /// colvars with multple items (e.g. multiple active CVCs) are duplicated
  std::vector<colvar *> colvars_smp;
  /// Indexes of the items to calculate for each colvar
  std::vector<int> colvars_smp_items;

  /// Array of named atom groups
  std::vector<atom_group *> named_atom_groups;
public:
  /// Register a named atom group into named_atom_groups
  inline void register_named_atom_group(atom_group * ag) {
    named_atom_groups.push_back(ag);
  }

  /// Array of collective variables
  std::vector<colvar *> *variables();

  /* TODO: implement named CVCs
  /// Array of named (reusable) collective variable components
  static std::vector<cvc *>     cvcs;
  /// Named cvcs register themselves at initialization time
  inline void register_cvc(cvc *p) {
    cvcs.push_back(p);
  }
  */

  /// Collective variables with the active flag on
  std::vector<colvar *> *variables_active();

  /// Collective variables to be calculated on different threads;
  /// colvars with multple items (e.g. multiple active CVCs) are duplicated
  std::vector<colvar *> *variables_active_smp();

  /// Indexes of the items to calculate for each colvar
  std::vector<int> *variables_active_smp_items();

  /// Array of collective variable biases
  std::vector<colvarbias *> biases;

  /// Energy of built-in and scripted biases, summed per time-step
  real total_bias_energy;

private:

  /// Array of active collective variable biases
  std::vector<colvarbias *> biases_active_;

public:

  /// Array of active collective variable biases
  std::vector<colvarbias *> *biases_active();

  /// \brief Whether debug output should be enabled (compile-time option)
  static inline bool debug()
  {
    return COLVARS_DEBUG;
  }

  /// \brief How many objects are configured yet?
  size_t size() const;

  /// \brief Constructor \param config_name Configuration file name
  /// \param restart_name (optional) Restart file name
  colvarmodule(colvarproxy *proxy);

  /// Destructor
  ~colvarmodule();

  /// Actual function called by the destructor
  int reset();

  /// Open a config file, load its contents, and pass it to config_string()
  int read_config_file(char const *config_file_name);

  /// \brief Parse a config string assuming it is a complete configuration
  /// (i.e. calling all parse functions)
  int read_config_string(std::string const &conf);

  /// \brief Parse a "clean" config string (no comments)
  int parse_config(std::string &conf);

  // Parse functions (setup internal data based on a string)

  /// Allow reading from Windows text files using using std::getline
  /// (which can still be used when the text is produced by Colvars itself)
  static std::istream & getline(std::istream &is, std::string &line);

  /// Parse the few module's global parameters
  int parse_global_params(std::string const &conf);

  /// Parse and initialize collective variables
  int parse_colvars(std::string const &conf);

  /// Parse and initialize collective variable biases
  int parse_biases(std::string const &conf);

  /// \brief Add new configuration during parsing (e.g. to implement
  /// back-compatibility); cannot be nested, i.e. conf should not contain
  /// anything that triggers another call
  int append_new_config(std::string const &conf);

private:

  /// Auto-generated configuration during parsing (e.g. to implement
  /// back-compatibility)
  std::string extra_conf;

  /// Parse and initialize collective variable biases of a specific type
  template <class bias_type>
  int parse_biases_type(std::string const &conf, char const *keyword);

  /// Test error condition and keyword parsing
  /// on error, delete new bias
  bool check_new_bias(std::string &conf, char const *key);

public:

  /// Return how many biases have this feature enabled
  static int num_biases_feature(int feature_id);

  /// Return how many biases are defined with this type
  static int num_biases_type(std::string const &type);

private:
  /// Useful wrapper to interrupt parsing if any error occurs
  int catch_input_errors(int result);

public:

  // "Setup" functions (change internal data based on related data
  // from the proxy that may change during program execution)
  // No additional parsing is done within these functions

  /// (Re)initialize internal data (currently used by LAMMPS)
  /// Also calls setup() member functions of colvars and biases
  int setup();

  /// (Re)initialize and (re)read the input state file calling read_restart()
  int setup_input();

  /// (Re)initialize the output trajectory and state file (does not write it yet)
  int setup_output();

#ifdef NAMD_VERSION
  typedef ofstream_namd ofstream;
#else
  typedef std::ofstream ofstream;
#endif

  /// Read the input restart file
  std::istream & read_restart(std::istream &is);
  /// Write the output restart file
  std::ostream & write_restart(std::ostream &os);

  /// Open a trajectory file if requested (and leave it open)
  int open_traj_file(std::string const &file_name);
  /// Close it
  int close_traj_file();
  /// Write in the trajectory file
  std::ostream & write_traj(std::ostream &os);
  /// Write explanatory labels in the trajectory file
  std::ostream & write_traj_label(std::ostream &os);

  /// Write all trajectory files
  int write_traj_files();
  /// Write all restart files
  int write_restart_files();
  /// Write all FINAL output files
  int write_output_files();
  /// Backup a file before writing it
  static int backup_file(char const *filename);

  /// Look up a bias by name; returns NULL if not found
  static colvarbias * bias_by_name(std::string const &name);

  /// Look up a colvar by name; returns NULL if not found
  static colvar * colvar_by_name(std::string const &name);

  /// Look up a named atom group by name; returns NULL if not found
  static atom_group * atom_group_by_name(std::string const &name);

  /// Load new configuration for the given bias -
  /// currently works for harmonic (force constant and/or centers)
  int change_configuration(std::string const &bias_name, std::string const &conf);

  /// Read a colvar value
  std::string read_colvar(std::string const &name);

  /// Calculate change in energy from using alt. config. for the given bias -
  /// currently works for harmonic (force constant and/or centers)
  real energy_difference(std::string const &bias_name, std::string const &conf);

  /// Give the total number of bins for a given bias.
  int bias_bin_num(std::string const &bias_name);
  /// Calculate the bin index for a given bias.
  int bias_current_bin(std::string const &bias_name);
  //// Give the count at a given bin index.
  int bias_bin_count(std::string const &bias_name, size_t bin_index);
  //// Share among replicas.
  int bias_share(std::string const &bias_name);

  /// Main worker function
  int calc();

  /// Calculate collective variables
  int calc_colvars();

  /// Calculate biases
  int calc_biases();

  /// Integrate bias and restraint forces, send colvar forces to atoms
  int update_colvar_forces();

  /// Perform analysis
  int analyze();

  /// \brief Read a collective variable trajectory (post-processing
  /// only, not called at runtime)
  int read_traj(char const *traj_filename,
                long        traj_read_begin,
                long        traj_read_end);

  /// Quick conversion of an object to a string
  template<typename T> static std::string to_str(T const &x,
                                                  size_t const &width = 0,
                                                  size_t const &prec = 0);
  /// Quick conversion of a vector of objects to a string
  template<typename T> static std::string to_str(std::vector<T> const &x,
                                                  size_t const &width = 0,
                                                  size_t const &prec = 0);

  /// Reduce the number of characters in a string
  static inline std::string wrap_string(std::string const &s,
                                         size_t const &nchars)
  {
    if (!s.size())
      return std::string(nchars, ' ');
    else
      return ( (s.size() <= size_t(nchars)) ?
               (s+std::string(nchars-s.size(), ' ')) :
               (std::string(s, 0, nchars)) );
  }

  /// Number of characters to represent a time step
  static size_t const it_width;
  /// Number of digits to represent a collective variables value(s)
  static size_t const cv_prec;
  /// Number of characters to represent a collective variables value(s)
  static size_t const cv_width;
  /// Number of digits to represent the collective variables energy
  static size_t const en_prec;
  /// Number of characters to represent the collective variables energy
  static size_t const en_width;
  /// Line separator in the log output
  static const char * const line_marker;


  // proxy functions

  /// \brief Value of the unit for atomic coordinates with respect to
  /// angstroms (used by some variables for hard-coded default values)
  static real unit_angstrom();

  /// \brief Boltmann constant
  static real boltzmann();

  /// \brief Temperature of the simulation (K)
  static real temperature();

  /// \brief Time step of MD integrator (fs)
  static real dt();

  /// Request calculation of total force from MD engine
  static void request_total_force();

  /// Print a message to the main log
  static void log(std::string const &message);

  /// Print a message to the main log and exit with error code
  static int fatal_error(std::string const &message);

  /// Print a message to the main log and set global error code
  static int error(std::string const &message, int code = COLVARS_ERROR);

  /// Print a message to the main log and exit normally
  static void exit(std::string const &message);

  // Replica exchange commands.
  static bool replica_enabled();
  static int replica_index();
  static int replica_num();
  static void replica_comm_barrier();
  static int replica_comm_recv(char* msg_data, int buf_len, int src_rep);
  static int replica_comm_send(char* msg_data, int msg_len, int dest_rep);

  /// \brief Get the distance between two atomic positions with pbcs handled
  /// correctly
  static rvector position_distance(atom_pos const &pos1,
                                    atom_pos const &pos2);


  /// \brief Get the square distance between two positions (with
  /// periodic boundary conditions handled transparently)
  ///
  /// Note: in the case of periodic boundary conditions, this provides
  /// an analytical square distance (while taking the square of
  /// position_distance() would produce leads to a cusp)
  static real position_dist2(atom_pos const &pos1,
                              atom_pos const &pos2);

  /// \brief Get the closest periodic image to a reference position
  /// \param pos The position to look for the closest periodic image
  /// \param ref_pos (optional) The reference position
  static void select_closest_image(atom_pos &pos,
                                    atom_pos const &ref_pos);

  /// \brief Perform select_closest_image() on a set of atomic positions
  ///
  /// After that, distance vectors can then be calculated directly,
  /// without using position_distance()
  static void select_closest_images(std::vector<atom_pos> &pos,
                                     atom_pos const &ref_pos);


  /// \brief Names of groups from a Gromacs .ndx file to be read at startup
  std::list<std::string> index_group_names;

  /// \brief Groups from a Gromacs .ndx file read at startup
  std::list<std::vector<int> > index_groups;

  /// \brief Read a Gromacs .ndx file
  int read_index_file(char const *filename);


  /// \brief Create atoms from a file \param filename name of the file
  /// (usually a PDB) \param atoms array of the atoms to be allocated
  /// \param pdb_field (optiona) if "filename" is a PDB file, use this
  /// field to determine which are the atoms to be set
  static int load_atoms(char const *filename,
                        atom_group &atoms,
                        std::string const &pdb_field,
                        double const pdb_field_value = 0.0);

  /// \brief Load the coordinates for a group of atoms from a file
  /// (PDB or XYZ)
  static int load_coords(char const *filename,
                         std::vector<atom_pos> &pos,
                         const std::vector<int> &indices,
                         std::string const &pdb_field,
                         double const pdb_field_value = 0.0);

  /// \brief Load the coordinates for a group of atoms from an
  /// XYZ file
  static int load_coords_xyz(char const *filename,
                              std::vector<atom_pos> &pos,
                              const std::vector<int> &indices);

  /// Frequency for collective variables trajectory output
  static size_t cv_traj_freq;

  /// \brief True if only analysis is performed and not a run
  static bool   b_analysis;

  /// Frequency for saving output restarts
  static size_t restart_out_freq;
  /// Output restart file name
  std::string   restart_out_name;

  /// Pseudo-random number with Gaussian distribution
  static real rand_gaussian(void);

protected:

  /// Configuration file
  std::ifstream config_s;

  /// Configuration file parser object
  colvarparse *parse;

  /// Name of the trajectory file
  std::string cv_traj_name;

  /// Collective variables output trajectory file
  colvarmodule::ofstream cv_traj_os;

  /// Appending to the existing trajectory file?
  bool cv_traj_append;

  /// Output restart file
  colvarmodule::ofstream restart_out_os;

private:

  /// Counter for the current depth in the object hierarchy (useg e.g. in output)
  size_t depth_s;

  /// Thread-specific depth
  std::vector<size_t> depth_v;

public:

  /// Get the current object depth in the hierarchy
  static size_t & depth();

  /// Increase the depth (number of indentations in the output)
  static void increase_depth();

  /// Decrease the depth (number of indentations in the output)
  static void decrease_depth();

  static inline bool scripted_forces()
  {
    return use_scripted_forces;
  }

  /// Use scripted colvars forces?
  static bool use_scripted_forces;

  /// Wait for all biases before calculating scripted forces?
  static bool scripting_after_biases;

  /// Calculate the energy and forces of scripted biases
  int calc_scripted_forces();

  /// \brief Pointer to the proxy object, used to retrieve atomic data
  /// from the hosting program; it is static in order to be accessible
  /// from static functions in the colvarmodule class
  static colvarproxy *proxy;

  /// \brief Accessor for the above
  static colvarmodule *main();

};


/// Shorthand for the frequently used type prefix
typedef colvarmodule cvm;


#include "colvartypes.h"


std::ostream & operator << (std::ostream &os, cvm::rvector const &v);
std::istream & operator >> (std::istream &is, cvm::rvector &v);


template<typename T> std::string cvm::to_str(T const &x,
                                              size_t const &width,
                                              size_t const &prec) {
  std::ostringstream os;
  if (width) os.width(width);
  if (prec) {
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(prec);
  }
  os << x;
  return os.str();
}

template<typename T> std::string cvm::to_str(std::vector<T> const &x,
                                              size_t const &width,
                                              size_t const &prec) {
  if (!x.size()) return std::string("");
  std::ostringstream os;
  if (prec) {
    os.setf(std::ios::scientific, std::ios::floatfield);
  }
  os << "{ ";
  if (width) os.width(width);
  if (prec) os.precision(prec);
  os << x[0];
  for (size_t i = 1; i < x.size(); i++) {
    os << ", ";
    if (width) os.width(width);
    if (prec) os.precision(prec);
    os << x[i];
  }
  os << " }";
  return os.str();
}


#include "colvarproxy.h"


inline cvm::real cvm::unit_angstrom()
{
  return proxy->unit_angstrom();
}

inline cvm::real cvm::boltzmann()
{
  return proxy->boltzmann();
}

inline cvm::real cvm::temperature()
{
  return proxy->temperature();
}

inline cvm::real cvm::dt()
{
  return proxy->dt();
}

// Replica exchange commands
inline bool cvm::replica_enabled() {
  return proxy->replica_enabled();
}
inline int cvm::replica_index() {
  return proxy->replica_index();
}
inline int cvm::replica_num() {
  return proxy->replica_num();
}
inline void cvm::replica_comm_barrier() {
  return proxy->replica_comm_barrier();
}
inline int cvm::replica_comm_recv(char* msg_data, int buf_len, int src_rep) {
  return proxy->replica_comm_recv(msg_data,buf_len,src_rep);
}
inline int cvm::replica_comm_send(char* msg_data, int msg_len, int dest_rep) {
  return proxy->replica_comm_send(msg_data,msg_len,dest_rep);
}


inline void cvm::request_total_force()
{
  proxy->request_total_force(true);
}

inline void cvm::select_closest_image(atom_pos &pos,
                                       atom_pos const &ref_pos)
{
  proxy->select_closest_image(pos, ref_pos);
}

inline void cvm::select_closest_images(std::vector<atom_pos> &pos,
                                        atom_pos const &ref_pos)
{
  proxy->select_closest_images(pos, ref_pos);
}

inline cvm::rvector cvm::position_distance(atom_pos const &pos1,
                                            atom_pos const &pos2)
{
  return proxy->position_distance(pos1, pos2);
}

inline cvm::real cvm::position_dist2(cvm::atom_pos const &pos1,
                                      cvm::atom_pos const &pos2)
{
  return proxy->position_dist2(pos1, pos2);
}

inline cvm::real cvm::rand_gaussian(void)
{
  return proxy->rand_gaussian();
}

#endif
