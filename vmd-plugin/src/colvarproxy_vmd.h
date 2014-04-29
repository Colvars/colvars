// -*- c++ -*-

#ifndef COLVARPROXY_VMD_H
#define COLVARPROXY_VMD_H

#include "DrawMolecule.h"
#include "Timestep.h"
#include "Inform.h"

#include "colvarmodule.h"
#include "colvarproxy.h"


/// \brief Communication between colvars and VMD (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_vmd : public colvarproxy {

protected:

  VMDApp *vmd;
  int vmdmolid;
  DrawMolecule *vmdmol;
  Inform msgColvars;

public:

  friend class cvm::atom;

  colvarproxy_vmd();
  ~colvarproxy_vmd();

  void update_conf();

  inline cvm::real unit_angstrom()
  {
    return 1.0;
  }

  cvm::real boltzmann()
  {
    return 0.001987191;
  }

  cvm::real temperature()
  {
    // TODO implement a user method to set the value of this
    return 300.0;
  }

  cvm::real dt()
  {
    // TODO implement a user method to set the value of this
    return 1.0;
  }

  cvm::real rand_gaussian()
  {
    return vmd_random_gaussian();
  }

  std::string input_prefix()
  {
    return std::string (vmdmol->molname());
  }

  std::string restart_output_prefix()
  {
    // note that this shouldn't be called while running VMD anyway
    return std::string (vmdmol->molname()) + std::string (".rst");
  }

  std::string output_prefix()
  {
    // note that this shouldn't be called while running VMD anyway
    return std::string (vmdmol->molname()) + std::string (".out");
  }

  size_t restart_frequency() {
    return 0;
  }

  void add_energy (cvm::real energy);

  /// nothing to do here
  inline void request_system_force (bool yesno) {}


  cvm::rvector position_distance (cvm::atom_pos const &pos1,
                                  cvm::atom_pos const &pos2);
  cvm::real position_dist2 (cvm::atom_pos const &pos1,
                            cvm::atom_pos const &pos2);

  void select_closest_image (cvm::atom_pos &pos,
                             cvm::atom_pos const &ref_pos);

  void log (std::string const &message);
  void fatal_error (std::string const &message);
  void exit (std::string const &message);


  void load_atoms (char const *filename,
                   std::vector<cvm::atom> &atoms,
                   std::string const pdb_field,
                   double const pdb_field_value = 0.0);

  void load_coords (char const *filename,
                    std::vector<cvm::atom_pos> &pos,
                    const std::vector<int> &indices,
                    std::string const pdb_field,
                    double const pdb_field_value = 0.0);

  // no need to reimplement backup_file()
};


inline cvm::rvector colvarproxy_vmd::position_distance (cvm::atom_pos const &pos1,
                                                        cvm::atom_pos const &pos2)
{
  // TODO: add in the proxy constructor a check for orthonormal PBCs
  Timestep *ts = mol->current();
  cvm::real const a = ts->a_length;
  cvm::real const b = ts->b_length;
  cvm::real const c = ts->c_length;
  cvm::rvector const diff = (pos2 - pos1);
  while (diff.x <= -0.5*a) diff.x += a;
  while (diff.y <= -0.5*b) diff.y += b;
  while (diff.z <= -0.5*c) diff.z += c;
  while (diff.x  >  0.5*a) diff.x -= a;
  while (diff.y  >  0.5*b) diff.y -= b;
  while (diff.z  >  0.5*c) diff.z -= c;
  return diff;
}


inline void colvarproxy_vmd::select_closest_image (cvm::atom_pos &pos,
                                                   cvm::atom_pos const &ref_pos)
{
  cvm::rvector const diff = position_distance (ref_pos, pos);
  pos = ref_pos + diff;
}


inline cvm::real colvarproxy_vmd::position_dist2 (cvm::atom_pos const &pos1,
                                                  cvm::atom_pos const &pos2)
{
  cvm::rvector const d = position_distance (pos1, pos2);
  return cvm::real (d.x*d.x + d.y*d.y + d.z*d.z);
}



#endif

