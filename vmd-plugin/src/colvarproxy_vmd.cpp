#include "DrawMolecule.h"
#include "Timestep.h"
#include "Residue.h"
#include "Inform.h"


#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "colvarproxy_vmd.h"

#if defined(VMDTKCON)
Inform msgColvars("colvars) ",    VMDCON_INFO);
#else
// XXX global instances of the Inform class
Inform msgColvars("colvars) ");
#endif


void colvarproxy_vmd::log (std::string const &message)
{
  std::istringstream is (message);
  std::string line;
  while (std::getline (is, line)) {
    msgColvars << line.c_str() << "\n";
  }
}


void colvarproxy_vmd::fatal_error (std::string const &message)
{
  cvm::log (message);
  if (!cvm::debug())
    cvm::log ("If this error message is unclear, "
              "try recompiling the colvars plugin with -DCOLVARS_DEBUG.\n");
  // TODO: return control to Tcl interpreter
}


void colvarproxy_vmd::exit (std::string const &message)
{
  cvm::log (message);
  // TODO: return control to Tcl interpreter
}



size_t colvarproxy_vmd::init_atom (int const &aid)
{
  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    if (colvars_atoms[i] == aid) {
      // this atom id was already recorded
      colvars_atoms_ncopies[i] += 1;
      return i;
    }
  }

  // allocate a new slot for this atom
  colvars_atoms_ncopies.push_back (1);
  colvars_atoms.push_back (aid);
  positions.push_back (cvm::rvector());
  total_forces.push_back (cvm::rvector());
  applied_forces.push_back (cvm::rvector());

  return (colvars_atoms.size()-1);
}


// atom member functions, VMD specific implementations

cvm::atom::atom (int const &atom_number)
{
  // VMD internal numbering starts from zero
  int const aid (atom_number-1);

  DrawMolecule *mol = ((colvarproxy_vmd *) cvm::proxy)->vmdmol;
  float *masses = mol->mass();

  if (cvm::debug())
    cvm::log ("Adding atom "+cvm::to_str (aid+1)+
              " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= mol->nAtoms) ) 
    cvm::fatal_error ("Error: invalid atom number specified, "+
                      cvm::to_str (atom_number)+"\n");


  this->index = ((colvarproxy_vmd *) cvm::proxy)->init_atom (aid);
  if (cvm::debug())
    cvm::log ("The index of this atom in the colvarproxy_vmd arrays is "+
              cvm::to_str (this->index)+".\n");
  this->id = aid;
  this->mass = masses[aid];
  this->reset_data();
}

int const lookup_atom_from_name (int resid, char *name, char *segname)
{ 

}



// In case of PSF structure, this function's argument "resid" is the non-unique identifier
// TODO: check that the default segment_id of non-PSF topologies is MAIN
cvm::atom::atom (cvm::residue_id const &resid,
                 std::string const     &atom_name,
                 std::string const     &segment_name)
{

  DrawMolecule *mol = ((colvarproxy_vmd *) cvm::proxy)->vmdmol;

  int aid = -1;
  for (int ir = 0; ir < mol->nResidues; ir++) {
    Residue *vmdres = mol.residue(ir);
    if (vmdres->resid == resid) {
      for (ia = 0; ia < vmdres->natoms; ia++) {
        int const resaid = vmdres->atoms[ia];
        std::string const sel_segname ((mol->segNames).name(mol->atom(resaid)->segnameindex));
        std::string const sel_atom_name ((mol->atomNames).name(mol->atom(resaid)->nameindex));
        if ( ((segment_name.size() == 0) || (segment_name == sel_segname)) &&
             (atom_name == sel_atom_name) ) {
          aid = resaid;
          break;
        }
      }
    }
    if (aid >= 0) break;
  }

  if (cvm::debug())
    cvm::log ("Adding atom \""+
              atom_name+"\" in residue "+
              cvm::to_str (residue)+
              " (index "+cvm::to_str (aid)+
              ") for collective variables calculation.\n");

  if (aid < 0) {
    cvm::fatal_error ("Error: could not find atom \""+
                      atom_name+"\" in residue "+
                      cvm::to_str (residue)+
                      ( (segment_name.size()) ?
                        (", segment \""+segment_name+"\"") :
                        ("") )+
                      "\n");
  }

  this->index = ((colvarproxy_vmd *) cvm::proxy)->init_atom (aid);
  if (cvm::debug())
    cvm::log ("The index of this atom in the colvarproxy_vmd arrays is "+
              cvm::to_str (this->index)+".\n");
  this->id = aid;
  this->mass = masses[aid];
  this->reset_data();
}


// copy constructor
cvm::atom::atom (cvm::atom const &a)
  : index (a.index), id (a.id), mass (a.mass)
{
  // increment the counter 
  colvarproxy_vmd *p = (colvarproxy_vmd *) cvm::proxy;
  p->colvars_atoms_ncopies[this->index] += 1;
}


cvm::atom::~atom() 
{
  if (this->index >= 0) {
    colvarproxy_vmd *p = (colvarproxy_vmd *) cvm::proxy;
    if (p->colvars_atoms_ncopies[this->index] > 0)
      p->colvars_atoms_ncopies[this->index] -= 1;
  }
}


void cvm::atom::read_position()
{
  // read the position directly from the current timestep's memory
  // Note: no prior update should be required (unlike NAMD with GlobalMaster)
  DrawMolecule *mol = ((colvarproxy_vmd *) cvm::proxy)->vmdmol;
  float *vmdpos = (mol->current())->vmdpos;
  this->pos = cvm::atom_pos (vmdpos[this->id*3+0],
                             vmdpos[this->id*3+1],
                             vmdpos[this->id*3+2]);
}
