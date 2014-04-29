#include "VMDApp.h"
#include "DrawMolecule.h"
#include "Timestep.h"
#include "Residue.h"
#include "Inform.h"
#include "utilities.h"

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "colvarproxy_vmd.h"


colvarproxy_vmd::colvarproxy_vmd (VMDApp *vmdapp)
  : vmd (vmdapp), vmdmolid (-1)
{
  first_timestep = true;
  system_force_requested = false;

  // same seed as in Measure.C
  vmd_srandom (38572111);

  vmdmolid = (vmd->molecule_top())->id();
  vmdmol = ((colvarproxy_vmd *) cvm::proxy)->vmd->moleculeList->mol_from_id (vmdmolid);

  update_conf();
}

void colvarproxy_vmd::update_conf()
{
  // TODO when implementing multiple instances
}


#if defined(VMDTKCON)
  Inform msgColvars("colvars) ",    VMDCON_INFO);
#else
  Inform msgColvars("colvars) ");
#endif

void colvarproxy_vmd::log (std::string const &message)
{
  std::istringstream is (message);
  std::string line;
  while (std::getline (is, line)) {
    msgColvars << line << "\n";
  }
}

void colvarproxy_vmd::fatal_error (std::string const &message)
{
  log (message);
  if (!cvm::debug())
    log ("If this error message is unclear, "
         "try recompiling the colvars plugin with -DCOLVARS_DEBUG.\n");
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
  // TODO: return control to Tcl interpreter
}

void colvarproxy_vmd::exit (std::string const &message)
{
  log (message);
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
  // TODO: return control to Tcl interpreter
}

void colvarproxy_vmd::add_energy (cvm::real energy)
{
  (vmdmol->current())->energy[TSE_RESTRAINT] += energy;
  (vmdmol->current())->energy[TSE_TOTAL] += energy;
}


e_pdb_field pdb_field_str2enum (std::string const &pdb_field_str)
{
  e_pdb_field pdb_field = e_pdb_none;

  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("O")) {
    pdb_field = e_pdb_occ;
  }

  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("B")) {
    pdb_field = e_pdb_beta;
  }

  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("X")) {
    pdb_field = e_pdb_x;
  }
  
  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("Y")) {
    pdb_field = e_pdb_y;
  }

  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("Z")) {
    pdb_field = e_pdb_z;
  }

  if (pdb_field == e_pdb_none) {
    cvm::fatal_error ("Error: unsupported PDB field, \""+
                      pdb_field_str+"\".\n");
  }

  return pdb_field;
}


void colvarproxy_vmd::load_coords (char const *pdb_filename,
                                   std::vector<cvm::atom_pos> &pos,
                                   const std::vector<int> &indices,
                                   std::string const pdb_field_str,
                                   double const pdb_field_value)
{
  if (pdb_field_str.size() == 0 && indices.size() == 0) {
    cvm::fatal_error ("Bug alert: either PDB field should be defined or list of "
                      "atom IDs should be available when loading atom coordinates!\n");
  }

  e_pdb_field pdb_field_index;
  bool const use_pdb_field = (pdb_field_str.size() > 0);
  if (use_pdb_field) {
    pdb_field_index = pdb_field_str2enum (pdb_field_str);
  }

  // next index to be looked up in PDB file (if list is supplied)
  std::vector<int>::const_iterator current_index = indices.begin();

  FileSpec *tmpspec = new FileSpec();
  int tmpmolid = vmd->molecule_load (-1, pdb_filename, "pdb", tmpspec);
  DrawMolecule *tmpmol = ((colvarproxy_vmd *) cvm::proxy)->vmd->moleculeList->mol_from_id (tmpmolid);
  delete tmpspec;
  vmd->molecule_make_top (vmdmolid);
  size_t const pdb_natoms = tmpmol->nAtoms;
  
  if (pos.size() != pdb_natoms) {

    bool const pos_allocated = (pos.size() > 0);

    size_t ipos = 0, ipdb = 0;
    for ( ; ipdb < pdb_natoms; ipdb++) {

      if (use_pdb_field) {
        // PDB field mode: skip atoms with wrong value in PDB field
        double atom_pdb_field_value = 0.0;

        switch (pdb_field_index) {
        case e_pdb_occ:
          atom_pdb_field_value = (tmpmol->occupancy())[ipdb];
          break;
        case e_pdb_beta:
          atom_pdb_field_value = (tmpmol->beta())[ipdb];
          break;
        case e_pdb_x:
          atom_pdb_field_value = (tmpmol->current()->pos)[ipdb*3];
          break;
        case e_pdb_y:
          atom_pdb_field_value = (tmpmol->current()->pos)[ipdb*3+1];
          break;
        case e_pdb_z:
          atom_pdb_field_value = (tmpmol->current()->pos)[ipdb*3+2];
          break;
        default:
          break;
        }

        if ( (pdb_field_value) &&
             (atom_pdb_field_value != pdb_field_value) ) {
          continue;
        } else if (atom_pdb_field_value == 0.0) {
          continue;
        }

      } else {
        // Atom ID mode: use predefined atom IDs from the atom group
        if (ipdb != *current_index) {
          // Skip atoms not in the list
          continue;
        } else {
          current_index++;
        }
      }
      
      if (!pos_allocated) {
        pos.push_back (cvm::atom_pos (0.0, 0.0, 0.0));
      } else if (ipos >= pos.size()) {
        cvm::fatal_error ("Error: the PDB file \""+
                          std::string (pdb_filename)+
                          "\" contains coordinates for "
                          "more atoms than needed.\n");
      }

      pos[ipos] = cvm::atom_pos ((tmpmol->current()->pos)[ipdb*3],
                                 (tmpmol->current()->pos)[ipdb*3+1],
                                 (tmpmol->current()->pos)[ipdb*3+2]);
      ipos++;
      if (!use_pdb_field && current_index == indices.end())
        break;
    }

    if ((ipos < pos.size()) || (current_index != indices.end()))
      cvm::fatal_error ("Error: the number of records in the PDB file \""+
                        std::string (pdb_filename)+
                        "\" does not appear to match either the total number of atoms,"+
                        " or the number of coordinates requested at this point ("+
                        cvm::to_str (pos.size())+").\n");

  } else {

    // when the PDB contains exactly the number of atoms of the array,
    // ignore the fields and just read coordinates
    for (size_t ia = 0; ia < pos.size(); ia++) {
      pos[ia] = cvm::atom_pos ((tmpmol->current()->pos)[ia*3],
                               (tmpmol->current()->pos)[ia*3+1],
                               (tmpmol->current()->pos)[ia*3+2]);
    }
  }

  vmd->molecule_delete (tmpmolid);
}



void colvarproxy_vmd::load_atoms (char const *pdb_filename,
                                  std::vector<cvm::atom> &atoms,
                                  std::string const pdb_field_str,
                                  double const pdb_field_value)
{
  if (pdb_field_str.size() == 0)
    cvm::fatal_error ("Error: must define which PDB field to use "
                      "in order to define atoms from a PDB file.\n");

  FileSpec *tmpspec = new FileSpec();
  int tmpmolid = vmd->molecule_load (-1, pdb_filename, "pdb", tmpspec);
  DrawMolecule *tmpmol = ((colvarproxy_vmd *) cvm::proxy)->vmd->moleculeList->mol_from_id (tmpmolid);
  delete tmpspec;
  vmd->molecule_make_top (vmdmolid);
  size_t const pdb_natoms = tmpmol->nAtoms;

  e_pdb_field pdb_field_index = pdb_field_str2enum (pdb_field_str);

  for (size_t ipdb = 0; ipdb < pdb_natoms; ipdb++) {

    double atom_pdb_field_value = 0.0;

    switch (pdb_field_index) {
    case e_pdb_occ:
      atom_pdb_field_value = (tmpmol->occupancy())[ipdb];
      break;
    case e_pdb_beta:
      atom_pdb_field_value = (tmpmol->beta())[ipdb];
      break;
    case e_pdb_x:
      atom_pdb_field_value = (tmpmol->current()->pos)[ipdb*3];
      break;
    case e_pdb_y:
      atom_pdb_field_value = (tmpmol->current()->pos)[ipdb*3+1];
      break;
    case e_pdb_z:
      atom_pdb_field_value = (tmpmol->current()->pos)[ipdb*3+2];
      break;
    default:
      break;
    }

    if ( (pdb_field_value) &&
         (atom_pdb_field_value != pdb_field_value) ) {
      continue;
    } else if (atom_pdb_field_value == 0.0) {
      continue;
    }
     
    atoms.push_back (cvm::atom (ipdb+1));
  }

  vmd->molecule_delete (tmpmolid);
}


// atom member functions, VMD specific implementations

cvm::atom::atom (int const &atom_number)
{
  // VMD internal numbering starts from zero
  int const aid (atom_number-1);

  float *masses = vmdmol->mass();

  if (cvm::debug())
    cvm::log ("Adding atom "+cvm::to_str (aid+1)+
              " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= vmdmol->nAtoms) ) 
    cvm::fatal_error ("Error: invalid atom number specified, "+
                      cvm::to_str (atom_number)+"\n");

  this->id = aid;
  this->mass = masses[aid];
  this->reset_data();
}


// In case of PSF structure, this function's argument "resid" is the non-unique identifier
// TODO: check that the default segment_id of non-PSF topologies is MAIN
cvm::atom::atom (cvm::residue_id const &resid,
                 std::string const     &atom_name,
                 std::string const     &segment_name)
{
  int aid = -1;
  for (int ir = 0; ir < vmdmol->nResidues; ir++) {
    Residue *vmdres = mol.residue(ir);
    if (vmdres->resid == resid) {
      for (ia = 0; ia < vmdres->natoms; ia++) {
        int const resaid = vmdres->atoms[ia];
        std::string const sel_segname ((vmdmol->segNames).name(vmdmol->atom(resaid)->segnameindex));
        std::string const sel_atom_name ((vmdmol->atomNames).name(vmdmol->atom(resaid)->nameindex));
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
  float *vmdpos = (vmdmol->current())->pos;
  this->pos = cvm::atom_pos (vmdpos[this->id*3+0],
                             vmdpos[this->id*3+1],
                             vmdpos[this->id*3+2]);
}
