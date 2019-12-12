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

protected:

  /// pointer to the VMD main object
  VMDApp *vmd;
  /// VMD molecule id being used (must be provided at construction)
  int vmdmolid;
  /// pointer to VMD molecule (derived from vmdmolid)
  DrawMolecule *vmdmol;
  /// current frame (returned by vmdmol->frame())
  long int vmdmol_frame;
  /// output object
  Inform msgColvars;

public:


  friend class cvm::atom;

  colvarproxy_vmd(Tcl_Interp *interp, VMDApp *vmd, int molid);
  ~colvarproxy_vmd();
  int request_deletion();

  int setup();

  int update_input();
  /// \brief Update mass, charge, etc
  int update_atomic_properties();

  inline cvm::real backend_angstrom_value()
  {
    return 1.0;
  }

  inline cvm::real boltzmann()
  {
    return 0.001987191;
  }

  inline cvm::real temperature()
  {
    // TODO define, document and implement a user method to set the value of this
    return 300.0;
  }

  inline cvm::real dt()
  {
    // TODO define, document and implement a user method to set the value of this
    return 1.0;
  }

  inline cvm::real rand_gaussian()
  {
    return vmd_random_gaussian();
  }

  /// Return molid of VMD molecule currently associated with Colvars
  inline int get_vmdmolid() { return vmdmolid; }

  inline int get_frame(long int &f)
  {
    f = vmdmol_frame;
    return COLVARS_OK;
  }

  int set_frame(long int f);

  std::string input_prefix_str;
  std::string input_prefix()
  {
    return input_prefix_str;
  }

  std::string restart_output_prefix()
  {
    // note: this shouldn't need to be called in VMD anyway
    return output_prefix_str;
  }

  std::string output_prefix_str;
  inline std::string output_prefix()
  {
    return output_prefix_str;
  }


#if defined(VMDTCL)
  void init_tcl_pointers();
#endif

  char const *script_obj_to_str(unsigned char *obj);
  std::vector<std::string> script_obj_to_str_vector(unsigned char *obj);

  void add_energy(cvm::real energy);

private:
  bool total_force_requested;
public:
  void request_total_force(bool yesno);

  std::string error_output;
  void log(std::string const &message);
  void error(std::string const &message);
  void fatal_error(std::string const &message);
  int set_unit_system(std::string const &units_in, bool check_only);

  // Callback functions
  int run_force_callback();
  int run_colvar_callback(std::string const &name,
                          std::vector<const colvarvalue *> const &cvcs,
                          colvarvalue &value);
  int run_colvar_gradient_callback(std::string const &name,
                                   std::vector<const colvarvalue *> const &cvc_values,
                                   std::vector<cvm::matrix2d<cvm::real> > &gradient);

  int load_atoms(char const *filename,
                 cvm::atom_group &atoms,
                 std::string const &pdb_field,
                 double const pdb_field_value = 0.0);

  int load_coords(char const *filename,
                  std::vector<cvm::atom_pos> &pos,
                  const std::vector<int> &indices,
                  std::string const &pdb_field,
                  double const pdb_field_value = 0.0);

  int init_atom(int atom_number);
  int check_atom_id(int atom_number);
  int init_atom(cvm::residue_id const &residue,
                std::string const     &atom_name,
                std::string const     &segment_id);
  int check_atom_id(cvm::residue_id const &residue,
                    std::string const     &atom_name,
                    std::string const     &segment_id);

};


#endif

