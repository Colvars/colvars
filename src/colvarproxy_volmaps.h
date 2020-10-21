// -*- c++ -*-

#ifndef COLVARPROXY_VOLMAPS_H
#define COLVARPROXY_VOLMAPS_H


/// \brief Container of grid-based objects
class colvarproxy_volmaps {

public:

  /// Contructor
  colvarproxy_volmaps();

  /// Destructor
  virtual ~colvarproxy_volmaps();

  /// Clear volumetric map data
  int reset();

  /// \brief Whether this implementation has capability to use volumetric maps
  virtual int volmaps_available();

  /// Create a slot for a volumetric map not requested yet
  int add_volmap_slot(int volmap_id);

  /// Request and prepare this volumetric map for use by Colvars
  /// \param volmap_id Numeric ID used by the MD engine
  /// \returns Index of the map in the colvarproxy arrays
  virtual int init_volmap_by_id(int volmap_id);

  /// Request and prepare this volumetric map for use by Colvars 
  /// \param volmap_name Name used by the MD engine
  /// \returns Index of the map in the colvarproxy arrays
  virtual int init_volmap_by_name(char const *volmap_name);

  /// Check that the given volmap ID is valid (return COLVARS_OK if it is)
  /// \param volmap_id Numeric ID used by the MD engine
  /// \returns Error code
  virtual int check_volmap_by_id(int volmap_id);

  /// Check that the given volmap name is valid (return COLVARS_OK if it is)
  /// \param volmap_name Name used by the MD engine
  /// \returns Error code
  virtual int check_volmap_by_name(char const *volmap_name);

  /// Request and prepare this volumetric map for use by Colvars
  int init_volmap_by_name(std::string const &volmap_name);

  /// Check that the given volmap name is valid (return COLVARS_OK if it is)
  int check_volmap_by_name(std::string const &volmap_name);

  /// \brief Used by the CVC destructors
  virtual void clear_volmap(int index);

  /// Get the numeric ID of the given volumetric map (for the MD program)
  virtual int get_volmap_id_from_name(char const *volmap_name);

  /// Get the numeric ID of the given volumetric map (for the MD program)
  inline int get_volmap_id(int index) const
  {
    return volmaps_ids[index];
  }

  /// Read the current value of the volumetric map
  inline cvm::real get_volmap_value(int index) const
  {
    return volmaps_values[index];
  }

  /// Request that this force is applied to the given volumetric map
  inline void apply_volmap_force(int index, cvm::real const &new_force)
  {
    volmaps_new_colvar_forces[index] += new_force;
  }

  /// Compute an atomic field (e.g. a colvar) weighted by a volumetric map
  /// \param volmap_id Numeric index of the map (no need to request it)
  /// \param atom_begin Iterator pointing to first atom
  /// \param atom_end Iterator pointing past the last atom
  /// \param value Pointer to location of total to increment
  /// \param atom_field Array of atomic field values (if NULL, 1 is used)
  virtual int compute_volmap(int volmap_id,
                             cvm::atom_iter atom_begin,
                             cvm::atom_iter atom_end,
                             cvm::real *value,
                             cvm::real *atom_field) const;

protected:

  /// \brief Array of numeric IDs of volumetric maps
  std::vector<int>          volmaps_ids;

  /// \brief Keep track of how many times each vol map is used by a
  /// separate colvar object
  std::vector<size_t>       volmaps_ncopies;

  /// \brief Current values of the vol maps
  std::vector<cvm::real>    volmaps_values;

  /// \brief Forces applied from colvars, to be communicated to the MD
  /// integrator
  std::vector<cvm::real>    volmaps_new_colvar_forces;
};


#endif
