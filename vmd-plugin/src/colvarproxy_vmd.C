#include <tcl.h>

#include "VMDApp.h"
#include "DrawMolecule.h"
#include "MoleculeList.h"
#include "Timestep.h"
#include "Residue.h"
#include "Inform.h"
#include "utilities.h"

#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "colvarproxy_vmd.h"


extern "C" {
  int tcl_colvars (ClientData nodata, Tcl_Interp *vmdtcl, int argc, const char *argv[]) {

    static colvarproxy_vmd *proxy = NULL;

    if (proxy != NULL) {

      if (proxy->script->args (argc, argv) == COLVARSCRIPT_ERROR) {
        return TCL_ERROR;
      } else {
        return TCL_OK;
      }
      
    } else {

      VMDApp *vmd = (VMDApp *) Tcl_GetAssocData (vmdtcl, "VMDApp", NULL);
      if (vmd == NULL) {
        Tcl_SetResult (vmdtcl, "Error: cannot find VMD main object.", TCL_STATIC);
        return TCL_ERROR;
      }

      if ( argc >= 3 ) {
        // require a molid to create the module
        if (!strcmp (argv[1], "molid")) {
          int molid = -1;
          if (!strcmp (argv[2], "top")) {
            molid = vmd->molecule_top();
          } else {
            Tcl_GetInt (vmdtcl, argv[2], &molid);
          }
          if (vmd->molecule_valid_id (molid)) {
            proxy = new colvarproxy_vmd (vmdtcl, vmd, molid);
            proxy->script = new colvarscript();
            return TCL_OK;
          }
        }
      }
    }

    Tcl_SetResult (vmdtcl, "usage: colvars molid <molecule id>", TCL_STATIC);
    return TCL_ERROR;
  }

  int Colvars_Init (Tcl_Interp *vmdtcl) {
    VMDApp *vmd = (VMDApp *) Tcl_GetAssocData (vmdtcl, "VMDApp", NULL);
    if (vmd == NULL) {
      Tcl_SetResult (vmdtcl, "Error: cannot find VMD main object.", TCL_STATIC);
      return TCL_ERROR;
    }
    Tcl_CreateCommand (vmdtcl, "colvars", tcl_colvars, (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL);
    Tcl_PkgProvide (vmdtcl, "colvars", COLVARS_VERSION);
    return TCL_OK;
  }
}


colvarproxy_vmd::colvarproxy_vmd (Tcl_Interp *vti, VMDApp *v, int molid)
  : vmdtcl (vti),
    vmd (v),
    vmdmolid (molid),
#if defined(VMDTKCON)
    msgColvars ("colvars: ",    VMDCON_INFO)
#else
    msgColvars ("colvars: ")
#endif
{
  colvars = NULL;

  // same seed as in Measure.C
  vmd_srandom (38572111);

  vmdmol = vmd->moleculeList->mol_from_id (vmdmolid);

  update_proxy_data();
}

void colvarproxy_vmd::update_proxy_data()
{
  // TODO when implementing multiple instances
}

void colvarproxy_vmd::log (std::string const &message)
{
  std::istringstream is (message);
  std::string line;
  while (std::getline (is, line)) {
    msgColvars << line.c_str() << sendmsg;
  }
}

void colvarproxy_vmd::fatal_error (std::string const &message)
{
  // TODO: return control to Tcl interpreter instead of exiting
  log (message);
  if (!cvm::debug())
    log ("If this error message is unclear, "
         "try recompiling the colvars plugin with -DCOLVARS_DEBUG.\n");
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
  vmd->VMDexit ("Collective variables error.\n", 1, 2);
}

void colvarproxy_vmd::exit (std::string const &message)
{
  // TODO: return control to Tcl interpreter
}

void colvarproxy_vmd::add_energy (cvm::real energy)
{
  (vmdmol->current())->energy[TSE_RESTRAINT] += energy;
  (vmdmol->current())->energy[TSE_TOTAL] += energy;
}


enum e_pdb_field {
  e_pdb_none,
  e_pdb_occ,
  e_pdb_beta,
  e_pdb_x,
  e_pdb_y,
  e_pdb_z,
  e_pdb_ntot
};

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
  DrawMolecule *tmpmol = vmd->moleculeList->mol_from_id (tmpmolid);
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
  DrawMolecule *tmpmol = vmd->moleculeList->mol_from_id (tmpmolid);
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

  DrawMolecule *vmdmol = ((colvarproxy_vmd *) cvm::proxy)->vmdmol;
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
  DrawMolecule *vmdmol = ((colvarproxy_vmd *) cvm::proxy)->vmdmol;
  float *masses = vmdmol->mass();
  
  int aid = -1;
  for (int ir = 0; ir < vmdmol->nResidues; ir++) {
    Residue *vmdres = vmdmol->residue (ir);
    if (vmdres->resid == resid) {
      for (int ia = 0; ia < vmdres->atoms.num(); ia++) {
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
              cvm::to_str (resid)+
              " (index "+cvm::to_str (aid)+
              ") for collective variables calculation.\n");

  if (aid < 0) {
    cvm::fatal_error ("Error: could not find atom \""+
                      atom_name+"\" in residue "+
                      cvm::to_str (resid)+
                      ( (segment_name.size()) ?
                        (", segment \""+segment_name+"\"") :
                        ("") )+
                      "\n");
  }

  this->id = aid;
  this->mass = masses[aid];
  this->reset_data();
}

