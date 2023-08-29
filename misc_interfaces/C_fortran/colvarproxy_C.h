// Test linking Colvars through a C interface

#include "colvarproxy.h"
#include <iostream>

/// \brief Communication between colvars and LAMMPS
/// (implementation of \link colvarproxy \endlink)
class colvarproxy_C : public colvarproxy {

public:
  colvarproxy_C();
  ~colvarproxy_C();

private:
  int set_unit_system(std::string const &units, bool check_only) { return 0.; }

  /// \brief Boltzmann constant
  cvm::real boltzmann() { return 0.; }

  /// \brief Target temperature of the simulation (K units)
  cvm::real temperature() { return 0.; }

  /// \brief Time step of the simulation (fs)
  cvm::real dt() { return 0.; }

  /// \brief Pseudo-random number with Gaussian distribution
  cvm::real rand_gaussian(void) { return 0.; }

  /// Pass restraint energy value for current timestep to MD engine
  void add_energy(cvm::real energy) {}

  /// numeric index (1-based)
  int init_atom(int atom_number) { return 0; }

  /// Check that this atom number is valid, but do not initialize the
  /// corresponding atom yet
  int check_atom_id(int atom_number) { return 0; }

  /// Print a message to the main log
  void log(std::string const &message) { std::cout << "colvars: " << message << std::endl; }

  /// Print a message to the main log and let the rest of the program handle the error
  void error(std::string const &message) { std::cout << "colvars: " << message << std::endl; }
};
