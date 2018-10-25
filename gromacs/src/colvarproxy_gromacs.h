/// -*- c++ -*-
// Jeff Comer's tests to see if he can link GROMACS and Colvars
#ifndef GMX_COLVARS_COLVARPROXY_GROMACS_H
#define GMX_COLVARS_COLVARPROXY_GROMACS_H

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/mdtypes/forceoutput.h"

#ifndef COLVARPROXY_VERSION
#define COLVARPROXY_VERSION "2018-10-05"
#endif

/// \brief Communication between colvars and Gromacs (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_gromacs : public colvarproxy {
public:
  // GROMACS structures.
  t_pbc gmx_pbc;
  t_mdatoms *gmx_atoms;
protected:
  cvm::real thermostat_temperature;
  cvm::real timestep;
  bool first_timestep;

  int outInd; // Index of the output files.
  bool colvars_restart;
  std::string config_file;
  size_t restart_frequency_s;
  int previous_gmx_step;
  bool total_force_requested;
  double bias_energy;

  // GROMACS random number generation.
  gmx::DefaultRandomEngine           rng;   // gromacs random number generator
  gmx::TabulatedNormalDistribution<> normal_distribution;

public:
  friend class cvm::atom;
  colvarproxy_gromacs();
  ~colvarproxy_gromacs();

  // Initialize colvars.
  void init(t_inputrec *gmx_inp, gmx_int64_t step);
  // Perform colvars computation, return bias energy.
  double calculate(gmx_int64_t step, const rvec *x, gmx::ForceWithVirial *force);
  void add_energy (cvm::real energy);

  // **************** SYSTEM-WIDE PHYSICAL QUANTITIES ****************
  cvm::real backend_angstrom_value();
  cvm::real boltzmann();
  cvm::real temperature();
  cvm::real dt();
  cvm::real rand_gaussian();
  // **************** SIMULATION PARAMETERS ****************
  void set_temper(double temper);
  size_t restart_frequency();
  std::string restart_output_prefix();
  std::string output_prefix();
  // **************** ACCESS ATOMIC DATA ****************
  void request_total_force (bool yesno);
  // **************** PERIODIC BOUNDARY CONDITIONS ****************
  cvm::rvector position_distance (cvm::atom_pos const &pos1,
				  cvm::atom_pos const &pos2);
  //cvm::real position_dist2 (cvm::atom_pos const &pos1,
  //                                    cvm::atom_pos const &pos2);
  void select_closest_image (cvm::atom_pos &pos,
                             cvm::atom_pos const &ref_pos);

  // **************** INPUT/OUTPUT ****************
  /// Print a message to the main log
  void log (std::string const &message);
  /// Print a message to the main log and let the rest of the program handle the error
  void error (std::string const &message);
  /// Print a message to the main log and exit with error code
  void fatal_error (std::string const &message);
  /// Print a message to the main log and exit normally
  void exit (std::string const &message);
  /// Request to set the units used internally by Colvars
  int set_unit_system(std::string const &units_in, bool check_only);
  int backup_file (char const *filename);
  /// Read atom identifiers from a file \param filename name of
  /// the file (usually a PDB) \param atoms array to which atoms read
  /// from "filename" will be appended \param pdb_field (optiona) if
  /// "filename" is a PDB file, use this field to determine which are
  /// the atoms to be set
  int load_atoms (char const *filename,
                           std::vector<cvm::atom> &atoms,
                           std::string const &pdb_field,
                           double const pdb_field_value = 0.0);
  /// Load the coordinates for a group of atoms from a file
  /// (usually a PDB); if "pos" is already allocated, the number of its
  /// elements must match the number of atoms in "filename"
  int load_coords (char const *filename,
                            std::vector<cvm::atom_pos> &pos,
                            const std::vector<int> &indices,
                            std::string const &pdb_field,
                            double const pdb_field_value = 0.0);

  int init_atom(int atom_number);

  int check_atom_id(int atom_number);
  void update_atom_properties(int index);
};

#endif
