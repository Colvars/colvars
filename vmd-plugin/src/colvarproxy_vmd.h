#ifndef COLVARPROXY_VMD_H
#define COLVARPROXY_VMD_H

#include "DrawMolecule.h"
#include "Timestep.h"

#include "colvarmodule.h"
#include "colvarproxy.h"


/// \brief Communication between colvars and VMD (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_vmd : public colvarproxy {

protected:

  DrawMolecule *vmdmol;

  std::string input_prefix_str, output_prefix_str, restart_output_prefix_str;
  size_t      restart_frequency_s;
  bool        first_timestep;
  bool        system_force_requested;

  std::vector<int>          colvars_atoms;
  std::vector<size_t>       colvars_atoms_ncopies;
  std::vector<cvm::rvector> positions;
  std::vector<cvm::rvector> total_forces;
  std::vector<cvm::rvector> applied_forces;

  size_t init_atom (in const &aid);

public:

  friend class cvm::atom;

  colvarproxy_vmd();
  ~colvarproxy_vmd();

  void log (std::string const &message);
  void fatal_error (std::string const &message);
  void exit (std::string const &message);

  inline cvm::real unit_angstrom()
  {
    return 1.0;
  }

  cvm::real boltzmann()
  {
    return 0.001987191;
  }

  inline std::string input_prefix()
  {
    return input_prefix_str;
  }

  cvm::rvector position_distance (cvm::atom_pos const &pos1,
                                  cvm::atom_pos const &pos2);
  cvm::real position_dist2 (cvm::atom_pos const &pos1,
                            cvm::atom_pos const &pos2);

  void select_closest_image (cvm::atom_pos &pos,
                             cvm::atom_pos const &ref_pos);


  void load_atoms (char const *filename,
                   std::vector<cvm::atom> &atoms,
                   std::string const pdb_field,
                   double const pdb_field_value = 0.0);

  void load_coords (char const *filename,
                    std::vector<cvm::atom_pos> &pos,
                    const std::vector<int> &indices,
                    std::string const pdb_field,
                    double const pdb_field_value = 0.0);

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


// Emacs
// Local Variables:
// mode: C++
// End:
