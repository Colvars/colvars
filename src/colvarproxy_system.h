// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_SYSTEM_H
#define COLVARPROXY_SYSTEM_H


/// Methods for accessing the simulation system (PBCs, integrator, etc)
class colvarproxy_system {

public:

  /// Constructor
  colvarproxy_system();

  /// Destructor
  virtual ~colvarproxy_system();

  /// \brief Name of the unit system used internally by Colvars (by default, that of the back-end).
  /// Supported depending on the back-end: real (A, kcal/mol), metal (A, eV), electron (Bohr, Hartree), gromacs (nm, kJ/mol)
  /// Note: calls to back-end PBC functions assume back-end length unit
  /// We use different unit from back-end in VMD bc using PBC functions from colvarproxy base class
  /// Colvars internal units are user specified, because the module exchanges info in unknown
  /// composite dimensions with user input, while it only exchanges quantities of known
  /// dimension with the back-end (length and forces)
  std::string units;

  /// \brief Request to set the units used internally by Colvars
  virtual int set_unit_system(std::string const &units, bool check_only);

  /// \brief Convert a length from Angstrom to internal
  inline cvm::real angstrom_to_internal(cvm::real l) const
  {
    return l * angstrom_value_;
  }

  /// \brief Convert a length from internal to Angstrom
  inline cvm::real internal_to_angstrom(cvm::real l) const
  {
    return l / angstrom_value_;
  }

  /// Boltzmann constant, with unit the same as energy / K
  inline cvm::real boltzmann() const
  {
    return boltzmann_;
  }

  /// Current target temperature of the simulation (K units)
  inline cvm::real target_temperature() const
  {
    return target_temperature_;
  }

  /// Set the current target temperature of the simulation (K units)
  virtual int set_target_temperature(cvm::real T);

  /// Time step of the simulation (fs units)
  inline double dt() const
  {
    return timestep_;
  }

  /// Set the current integration timestep of the simulation (fs units)
  virtual int set_integration_timestep(cvm::real dt);

  /// Time step of the simulation (fs units)
  inline int time_step_factor() const
  {
    return time_step_factor_;
  }

  /// Set the current integration timestep of the simulation (fs units)
  virtual int set_time_step_factor(int fact);

  /// \brief Pseudo-random number with Gaussian distribution
  virtual cvm::real rand_gaussian(void);

  /// Pass restraint energy value for current timestep to MD engine
  virtual void add_energy(cvm::real energy);

  /// Use the PBC functions from the Colvars library (as opposed to MD engine)
  inline bool & use_internal_pbc() { return use_internal_pbc_; }

  /// Get the PBC-aware distance vector between two positions (using the MD engine's convention)
  virtual cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                         cvm::atom_pos const &pos2) const;

  /// Inline version of position_distance()
  cvm::rvector position_distance_internal(cvm::atom_pos const &pos1,
                                          cvm::atom_pos const &pos2) const;

  /// Kernel used by position_distance_internal()
  /// @param pos1 First position
  /// @param pos2 Second position
  /// @param a Unit cell vector
  /// @param b Unit cell vector
  /// @param c Unit cell vector
  /// @param a_r Reciprocal cell vector
  /// @param b_r Reciprocal cell vector
  /// @param c_r Reciprocal cell vector
  /// @param a_p Is this dimension periodic?
  /// @param b_p Is this dimension periodic?
  /// @param c_p Is this dimension periodic?
  /// @return (pos2 - pos1) without periodicity, minimum-image difference otherwise
  static cvm::rvector position_distance_kernel(cvm::atom_pos const &pos1,
                                               cvm::atom_pos const &pos2,
                                               cvm::rvector const &a,
                                               cvm::rvector const &b,
                                               cvm::rvector const &c,
                                               cvm::rvector const &a_r,
                                               cvm::rvector const &b_r,
                                               cvm::rvector const &c_r,
                                               bool a_p = false,
                                               bool b_p = false,
                                               bool c_p = false);

  /// Recompute PBC reciprocal lattice (assumes XYZ periodicity)
  void update_pbc_lattice();

  /// Set the lattice vectors to zero
  void reset_pbc_lattice();

  /// \brief Tell the proxy whether total forces are needed (they may not
  /// always be available)
  virtual void request_total_force(bool yesno);

  /// Are total forces being used?
  virtual bool total_forces_enabled() const;

  /// Are total forces from the current step available?
  /// in which case they are really system forces
  virtual bool total_forces_same_step() const;

  /// Get the molecule ID when called in VMD; raise error otherwise
  /// \param molid Set this argument equal to the current VMD molid
  virtual int get_molid(int &molid);

  /// Get value of alchemical lambda parameter from back-end (if available)
  virtual int get_alch_lambda(cvm::real* lambda);

  /// Set value of alchemical lambda parameter to be sent to back-end at end of timestep
  void set_alch_lambda(cvm::real lambda);

  /// Send cached value of alchemical lambda parameter to back-end (if available)
  virtual int send_alch_lambda();

  /// Request energy computation every freq steps (necessary for NAMD3, not all back-ends)
  virtual int request_alch_energy_freq(int const /* freq */) {
    return COLVARS_OK;
  }

  /// Get energy derivative with respect to lambda (if available)
  virtual int get_dE_dlambda(cvm::real* dE_dlambda);

  /// Apply a scalar force on dE_dlambda (back-end distributes it onto atoms)
  virtual int apply_force_dE_dlambda(cvm::real* force);

  /// Get energy second derivative with respect to lambda (if available)
  virtual int get_d2E_dlambda2(cvm::real* d2E_dlambda2);

  /// Force to be applied onto alch. lambda, propagated from biasing forces on dE_dlambda
  cvm::real indirect_lambda_biasing_force;

  /// Get weight factor from accelMD
  virtual cvm::real get_accelMD_factor() const {
    cvm::error("Error: accessing the reweighting factor of accelerated MD  "
               "is not yet implemented in the MD engine.\n",
               COLVARS_NOT_IMPLEMENTED);
    return 1.0;
  }
  virtual bool accelMD_enabled() const {
    return false;
  }

protected:

  /// Next value of lambda to be sent to back-end
  cvm::real cached_alch_lambda;

  /// Whether lambda has been set and needs to be updated in backend
  bool cached_alch_lambda_changed;

  /// Boltzmann constant in internal Colvars units
  cvm::real boltzmann_;

  /// Most up to date target temperature (K units); default to 0.0 if undefined
  cvm::real target_temperature_;

  /// Current integration timestep (engine units); default to 1.0 if undefined
  double timestep_;

  /// Current timestep multiplier, if Colvars is only called once every n MD timesteps
  int time_step_factor_ = 1;

  /// \brief Value of 1 Angstrom in the internal (front-end) Colvars unit for atomic coordinates
  /// * defaults to 0 in the base class; derived proxy classes must set it
  /// * in VMD proxy, can only be changed when no variables are defined
  /// as user-defined values in composite units must be compatible with that system
  cvm::real angstrom_value_;

  /// \brief Value of 1 kcal/mol in the internal Colvars unit for energy
  cvm::real kcal_mol_value_;

  /// Whether the total forces have been requested
  bool total_force_requested;

  /// Use the PBC functions from the Colvars library (as opposed to MD engine)
  bool use_internal_pbc_ = false;

  /// Type of boundary conditions defined for the current computation
  ///
  /// Orthogonal and triclinic cells are made available to objects.
  /// For any other conditions (mixed periodicity, triclinic cells in LAMMPS)
  /// minimum-image distances are computed by the host engine by default
  enum Boundaries_type {
    boundaries_non_periodic,
    boundaries_pbc_ortho,
    boundaries_pbc_triclinic,
    boundaries_unsupported
  };

  /// Type of boundary conditions
  Boundaries_type boundaries_type;

  /// Bravais lattice vectors
  cvm::rvector unit_cell_x, unit_cell_y, unit_cell_z;

  /// Reciprocal lattice vectors
  cvm::rvector reciprocal_cell_x, reciprocal_cell_y, reciprocal_cell_z;
};


inline cvm::rvector colvarproxy_system::position_distance_internal(cvm::atom_pos const &pos1,
                                                                   cvm::atom_pos const &pos2) const
{
  if (boundaries_type == boundaries_unsupported) {
    cvm::error("Error: unsupported boundary conditions.\n", COLVARS_INPUT_ERROR);
    return cvm::rvector({0.0, 0.0, 0.0});
  }

  if (boundaries_type == boundaries_non_periodic) {
    return pos2 - pos1;
  }

  // Periodicity flags are hard-coded, because this is the only case supported so far other than
  // the two above
  return position_distance_kernel(pos1, pos2, unit_cell_x, unit_cell_y, unit_cell_z,
                                  reciprocal_cell_x, reciprocal_cell_y, reciprocal_cell_z,
                                  true, true, true);
}


inline cvm::rvector colvarproxy_system::position_distance_kernel(cvm::atom_pos const &pos1,
                                                                 cvm::atom_pos const &pos2,
                                                                 cvm::rvector const &a,
                                                                 cvm::rvector const &b,
                                                                 cvm::rvector const &c,
                                                                 cvm::rvector const &a_r,
                                                                 cvm::rvector const &b_r,
                                                                 cvm::rvector const &c_r,
                                                                 bool a_p,
                                                                 bool b_p,
                                                                 bool c_p)
{
  cvm::rvector diff = (pos2 - pos1);

  cvm::real const x_shift = std::floor(a_r * diff + 0.5);
  cvm::real const y_shift = std::floor(b_r * diff + 0.5);
  cvm::real const z_shift = std::floor(c_r * diff + 0.5);

  if (a_p) {
    diff.x -= x_shift * a.x + y_shift * b.x + z_shift * c.x;
  }

  if (b_p) {
    diff.y -= x_shift * a.y + y_shift * b.y + z_shift * c.y;
  }

  if (c_p) {
    diff.z -= x_shift * a.z + y_shift * b.z + z_shift * c.z;
  }

  return diff;
}


#endif
