/// -*- c++ -*-

#ifndef COLVARPROXY_VMD_H
#define COLVARPROXY_VMD_H

#include <tcl.h>

#include "DrawMolecule.h"
#include "Timestep.h"
#include "Inform.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvartypes.h"
#include "colvaratoms.h"

#ifndef COLVARPROXY_VERSION
#define COLVARPROXY_VERSION "2015-11-02"
#endif


int tcl_colvars(ClientData clientdata, Tcl_Interp *interp, int argc, const char *argv[]);


/// \brief Communication between colvars and VMD (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_vmd : public colvarproxy {

protected:

  /// pointer to the VMD Tcl interpreter
  Tcl_Interp *interp;
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

  int setup();

  int update_input();
  /// \brief Update mass, charge, etc
  int update_atomic_properties();

  inline cvm::real unit_angstrom()
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

  void add_energy(cvm::real energy);

private:
  bool total_force_requested;
public:
  void request_total_force(bool yesno);

  cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                  cvm::atom_pos const &pos2);
  cvm::real position_dist2(cvm::atom_pos const &pos1,
                            cvm::atom_pos const &pos2);

  void select_closest_image(cvm::atom_pos &pos,
                             cvm::atom_pos const &ref_pos);

  std::string error_output;
  void log(std::string const &message);
  void error(std::string const &message);
  void fatal_error(std::string const &message);
  void exit(std::string const &message);

  // Callback functions
  int run_force_callback();
  int run_colvar_callback(std::string const &name,
                          std::vector<const colvarvalue *> const &cvcs,
                          colvarvalue &value);
  int run_colvar_gradient_callback(std::string const &name,
                                   std::vector<const colvarvalue *> const &cvcs,
                                   std::vector<colvarvalue> &gradient);

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



inline cvm::rvector colvarproxy_vmd::position_distance(cvm::atom_pos const &pos1,
                                                       cvm::atom_pos const &pos2)
{
  // TODO: add in the proxy constructor a check for orthonormal PBCs
  Timestep *ts = vmdmol->get_frame(vmdmol_frame);
  cvm::real const a = ts->a_length;
  cvm::real const b = ts->b_length;
  cvm::real const c = ts->c_length;
  cvm::rvector diff = (pos2 - pos1);
  if (a*a > 1.0e-12) {
    while (diff.x <= -0.5*a) diff.x += a;
    while (diff.x  >  0.5*a) diff.x -= a;
  }
  if (b*b > 1.0e-12) {
    while (diff.y <= -0.5*b) diff.y += b;
    while (diff.y  >  0.5*b) diff.y -= b;
  }
  if (c*c > 1.0e-12) {
    while (diff.z <= -0.5*c) diff.z += c;
    while (diff.z  >  0.5*c) diff.z -= c;
  }
  return diff;
}


inline void colvarproxy_vmd::select_closest_image(cvm::atom_pos &pos,
                                                  cvm::atom_pos const &ref_pos)
{
  cvm::rvector const diff = position_distance(ref_pos, pos);
  pos = ref_pos + diff;
}


inline cvm::real colvarproxy_vmd::position_dist2(cvm::atom_pos const &pos1,
                                                 cvm::atom_pos const &pos2)
{
  cvm::rvector const d = position_distance(pos1, pos2);
  return cvm::real(d.x*d.x + d.y*d.y + d.z*d.z);
}



#endif

