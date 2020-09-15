// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_VMD_H
#define COLVARPROXY_VMD_H

#include "colvarproxy_vmd_version.h"

#include <tcl.h>

#include "DrawMolecule.h"
#include "Timestep.h"
#include "Inform.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvartypes.h"
#include "colvaratoms.h"


int tcl_colvars(ClientData clientData, Tcl_Interp *interp,
                int objc, Tcl_Obj *const objv[]);

/// \brief Communication between colvars and VMD (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_vmd : public colvarproxy {

public:

  /// Constructor
  /// \param interl Pointer to Tcl interpreter
  /// \param vmd Pointer to VMDApp object
  /// \param molid Molecule ID (>= 0)
  colvarproxy_vmd(Tcl_Interp *interp, VMDApp *vmd, int molid);

  /// \brief Update mass, charge, etc
  int update_atomic_properties();

  virtual ~colvarproxy_vmd();

  virtual int request_deletion();

  virtual int setup();

  virtual int update_input();

  virtual cvm::real backend_angstrom_value();

  virtual cvm::real boltzmann();

  virtual cvm::real temperature();

  virtual cvm::real dt();

  virtual cvm::real rand_gaussian();

  virtual int get_molid(int &molid);

  virtual int get_frame(long int &f);

  virtual int set_frame(long int f);

  virtual void init_tcl_pointers();

  virtual void add_energy(cvm::real energy);

  virtual void request_total_force(bool yesno);

  virtual void log(std::string const &message);

  virtual void error(std::string const &message);

  virtual int set_unit_system(std::string const &units_in, bool check_only);

  virtual int run_force_callback();

  virtual int run_colvar_callback(std::string const &name,
                                  std::vector<const colvarvalue *> const &cvcs,
                                  colvarvalue &value);

  virtual int run_colvar_gradient_callback(std::string const &name,
                                           std::vector<const colvarvalue *> const &cvc_values,
                                           std::vector<cvm::matrix2d<cvm::real> > &gradient);

  virtual int load_atoms(char const *filename,
                         cvm::atom_group &atoms,
                         std::string const &pdb_field,
                         double const pdb_field_value = 0.0);

  virtual int load_coords(char const *filename,
                          std::vector<cvm::atom_pos> &pos,
                          const std::vector<int> &indices,
                          std::string const &pdb_field,
                          double const pdb_field_value = 0.0);

  virtual int init_atom(int atom_number);

  virtual int check_atom_id(int atom_number);

  virtual int init_atom(cvm::residue_id const &residue,
                        std::string const     &atom_name,
                        std::string const     &segment_id);

  virtual int check_atom_id(cvm::residue_id const &residue,
                            std::string const     &atom_name,
                            std::string const     &segment_id);

  virtual int init_volmap_by_id(int volmap_id);

  virtual int check_volmap_by_id(int volmap_id);

  virtual void clear_volmap(int index);

  virtual int compute_volmap(int flags,
                             int volmap_id,
                             cvm::atom_iter atom_begin,
                             cvm::atom_iter atom_end,
                             cvm::real *value,
                             cvm::real *atom_field);

  template<int flags>
  void compute_voldata(VolumetricData const *voldata,
                       cvm::atom_iter atom_begin,
                       cvm::atom_iter atom_end,
                       cvm::real *value,
                       cvm::real *atom_field);

protected:

  /// pointer to the VMD main object
  VMDApp *vmd;

  /// VMD molecule ID being used (must be provided at construction)
  int vmdmolid;

  /// pointer to VMD molecule (derived from vmdmolid)
  DrawMolecule *vmdmol;

  /// current frame (returned by vmdmol->frame())
  long int vmdmol_frame;

  /// output object
  Inform msgColvars;
};


#endif
