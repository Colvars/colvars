/// -*- c++ -*-

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
#include "colvarscript.h"
#include "colvarproxy.h"
#include "colvarproxy_vmd.h"


int tcl_colvars(ClientData clientdata, Tcl_Interp *interp, int argc, const char *argv[]) {

  static colvarproxy_vmd *proxy = NULL;
  int retval;

  if (proxy != NULL) {

    if ( argc >= 3 ) {
      if (!strcmp(argv[1], "molid")) {
        Tcl_SetResult(interp, (char *) (std::string("Colvars module already created: type \"cv\" for a list of arguments.").c_str()), TCL_STATIC);
        return TCL_ERROR;
      }
    }

    // Clear non-fatal errors from previous commands
    cvm::clear_error();

    retval = proxy->script->run(argc, argv);
    Tcl_SetResult(interp, (char *) proxy->script->result.c_str(), TCL_STATIC);

    if (cvm::get_error() & DELETE_COLVARS) {
      delete proxy;
      proxy = NULL;
      return TCL_OK;
    }

    if (cvm::get_error() & FATAL_ERROR) {
      // Fatal error: clean up cvm object and proxy
      delete proxy;
      proxy = NULL;
      return TCL_ERROR;
    }

    if (retval == COLVARSCRIPT_OK && !cvm::get_error()) {
      return TCL_OK;
    } else {
      return TCL_ERROR;
    }

  } else {

    VMDApp *vmd = (VMDApp *) clientdata;
    if (vmd == NULL) {
      Tcl_SetResult(interp, (char *) (std::string("Error: cannot find VMD main object.").c_str()), TCL_STATIC);
      return TCL_ERROR;
    }

    if ( argc >= 3 ) {
      // require a molid to create the module
      if (!strcmp(argv[1], "molid")) {
        int molid = -1;
        if (!strcmp(argv[2], "top")) {
          molid = vmd->molecule_top();
        } else {
          Tcl_GetInt(interp, argv[2], &molid);
        }
        if (vmd->molecule_valid_id(molid)) {
          proxy = new colvarproxy_vmd(interp, vmd, molid);
          return TCL_OK;
        } else {
          Tcl_SetResult(interp, (char *) (std::string("Error: molecule not found.").c_str()), TCL_STATIC);
          return TCL_ERROR;
        }
      }
    }
  }

  Tcl_SetResult(interp, (char *) (std::string("First, setup the colvars module with: cv molid <molecule id>").c_str()), TCL_STATIC);
  return TCL_ERROR;
}


colvarproxy_vmd::colvarproxy_vmd(Tcl_Interp *vti, VMDApp *v, int molid)
  : interp(vti),
    vmd(v),
    vmdmolid(molid),
#if defined(VMDTKCON)
    msgColvars("colvars: ",    VMDCON_INFO),
#else
    msgColvars("colvars: "),
#endif
    input_prefix_str(""), output_prefix_str("")
{
  // The module is only allocated here: it will be configured
  // through the "configfile" and "configstring" commands of colvarscript.
  colvars = new colvarmodule(this);
  cvm::log("Using VMD interface, version "+
            cvm::to_str(COLVARPROXY_VERSION)+".\n");

  colvars->cv_traj_freq = 0;
  colvars->restart_out_freq = 0;
  cvm::rotation::monitor_crossings = false;

  colvars->setup_input();
  colvars->setup_output();

  script = new colvarscript(this);
  script->proxy_error = COLVARSCRIPT_OK;

  // Do we have scripts?
  // For now colvars depend on Tcl, but this may not always be the case in the future
#if defined(VMDTCL)
  have_scripts = true;

  // User-scripted forces are not really useful in VMD, but we accept them
  // for compatibility with NAMD scripts
  if (Tcl_FindCommand(interp, "calc_colvar_forces", NULL, 0) == NULL) {
    force_script_defined = false;
  } else {
    force_script_defined = true;
  }
#else
  have_scripts = false;
#endif

  this->setup();
}


colvarproxy_vmd::~colvarproxy_vmd()
{
  if (script != NULL) {
    delete script;
    script = NULL;
  }
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
}


void colvarproxy_vmd::setup()
{
  vmdmol = vmd->moleculeList->mol_from_id(vmdmolid);
  if (vmdmol) {
    vmdmol_frame = vmdmol->frame();
  } else {
    fatal_error("Error: cannot find the molecule requested("+cvm::to_str(vmdmolid)+").\n");
  }
  if (colvars) {
    colvars->setup();
  }

  // same seed as in Measure.C
  vmd_srandom(38572111);
}


void colvarproxy_vmd::log(std::string const &message)
{
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line)) {
    msgColvars << line.c_str() << sendmsg;
  }
}

void colvarproxy_vmd::error(std::string const &message)
{
  log(message);
}

void colvarproxy_vmd::fatal_error(std::string const &message)
{
  log(message);
  if (!cvm::debug())
    log("If this error message is unclear, "
         "try recompiling VMD with -DCOLVARS_DEBUG.\n");
}

void colvarproxy_vmd::exit(std::string const &message)
{
  // Ultimately, this should never be called
  vmd->VMDexit("Collective variables module requested VMD shutdown.\n", 0, 2);
}


// Callback functions

#ifdef VMDTCL
int colvarproxy_vmd::run_force_callback() {
  std::string cmd = std::string("calc_colvar_forces ")
    + cvm::to_str(cvm::step_absolute());
  int err = Tcl_Eval(interp, cmd.c_str());
  if (err != TCL_OK) {
    cvm::log(std::string("Error while executing calc_colvar_forces:\n"));
    cvm::error(Tcl_GetStringResult(interp));
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}

int colvarproxy_vmd::run_colvar_callback(std::string const &name,
                      std::vector<const colvarvalue *> const &cvc_values,
                      colvarvalue &value)
{
  size_t i;
  std::string cmd = std::string("calc_") + name;
  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") +  (*(cvc_values[i])).to_simple_string() + std::string("}");
  }
  int err = Tcl_Eval(interp, cmd.c_str());
  const char *result = Tcl_GetStringResult(interp);
  if (err != TCL_OK) {
    cvm::log(std::string("Error while executing ")
              + cmd + std::string(":\n"));
    cvm::error(result);
    return COLVARS_ERROR;
  }
  std::istringstream is(result);
  if (value.from_simple_string(is.str()) != COLVARS_OK) {
    cvm::log("Error parsing colvar value from script:");
    cvm::error(result);
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}

int colvarproxy_vmd::run_colvar_gradient_callback(std::string const &name,
                               std::vector<const colvarvalue *> const &cvc_values,
                               std::vector<colvarvalue> &gradient)
{
  size_t i;
  std::string cmd = std::string("calc_") + name + "_gradient";
  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") +  (*(cvc_values[i])).to_simple_string() + std::string("}");
  }
  int err = Tcl_Eval(interp, cmd.c_str());
  if (err != TCL_OK) {
    cvm::log(std::string("Error while executing ")
              + cmd + std::string(":\n"));
    cvm::error(Tcl_GetStringResult(interp));
    return COLVARS_ERROR;
  }
  Tcl_Obj **list;
  int n;
  Tcl_ListObjGetElements(interp, Tcl_GetObjResult(interp),
                         &n, &list);
  if (n != int(gradient.size())) {
    cvm::error("Error parsing list of gradient values from script");
    return COLVARS_ERROR;
  }
  for (i = 0; i < gradient.size(); i++) {
    std::istringstream is(Tcl_GetString(list[i]));
    gradient[i].type(*(cvc_values[i]));
    gradient[i].is_derivative();
    if (gradient[i].from_simple_string(is.str()) != COLVARS_OK) {
      cvm::error("Error parsing gradient value from script");
      return COLVARS_ERROR;
    }
  }
  return (err == TCL_OK) ? COLVARS_OK : COLVARS_ERROR;
}
#endif

void colvarproxy_vmd::add_energy(cvm::real energy)
{
  // (vmdmol->get_frame (this->frame()))->energy[TSE_RESTRAINT] += energy;
  // (vmdmol->get_frame (this->frame()))->energy[TSE_TOTAL] += energy;
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

e_pdb_field pdb_field_str2enum(std::string const &pdb_field_str)
{
  e_pdb_field pdb_field = e_pdb_none;

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("O")) {
    pdb_field = e_pdb_occ;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("B")) {
    pdb_field = e_pdb_beta;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("X")) {
    pdb_field = e_pdb_x;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("Y")) {
    pdb_field = e_pdb_y;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("Z")) {
    pdb_field = e_pdb_z;
  }

  if (pdb_field == e_pdb_none) {
    cvm::fatal_error("Error: unsupported PDB field, \""+
                      pdb_field_str+"\".\n");
  }

  return pdb_field;
}


int colvarproxy_vmd::load_coords(char const *pdb_filename,
                                   std::vector<cvm::atom_pos> &pos,
                                   const std::vector<int> &indices,
                                   std::string const &pdb_field_str,
                                   double const pdb_field_value)
{
  if (pdb_field_str.size() == 0 && indices.size() == 0) {
    cvm::error("Bug alert: either PDB field should be defined or list of "
                "atom IDs should be available when loading atom coordinates!\n",
                BUG_ERROR);
    return COLVARS_ERROR;
  }

  e_pdb_field pdb_field_index;
  bool const use_pdb_field = (pdb_field_str.size() > 0);
  if (use_pdb_field) {
    pdb_field_index = pdb_field_str2enum(pdb_field_str);
  }

  // next index to be looked up in PDB file (if list is supplied)
  std::vector<int>::const_iterator current_index = indices.begin();

  FileSpec *tmpspec = new FileSpec();
  tmpspec->autobonds = 0;
  int tmpmolid = vmd->molecule_load(-1, pdb_filename, "pdb", tmpspec);
  delete tmpspec;
  if (tmpmolid < 0) {
    cvm::error("Error: VMD could not read file \""+std::string(pdb_filename)+"\".\n",
                FILE_ERROR);
    return COLVARS_ERROR;
  }
  DrawMolecule *tmpmol = vmd->moleculeList->mol_from_id(tmpmolid);

  vmd->molecule_make_top(vmdmolid);
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
          atom_pdb_field_value = (tmpmol->get_frame(0)->pos)[ipdb*3];
          break;
        case e_pdb_y:
          atom_pdb_field_value = (tmpmol->get_frame(0)->pos)[ipdb*3+1];
          break;
        case e_pdb_z:
          atom_pdb_field_value = (tmpmol->get_frame(0)->pos)[ipdb*3+2];
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
        if (((int)ipdb) != *current_index) {
          // Skip atoms not in the list
          continue;
        } else {
          current_index++;
        }
      }

      if (!pos_allocated) {
        pos.push_back(cvm::atom_pos(0.0, 0.0, 0.0));
      } else if (ipos >= pos.size()) {
        cvm::error("Error: the PDB file \""+
                    std::string(pdb_filename)+
                    "\" contains coordinates for "
                    "more atoms than needed.\n", INPUT_ERROR);
        return COLVARS_ERROR;
      }

      pos[ipos] = cvm::atom_pos((tmpmol->get_frame(0)->pos)[ipdb*3],
                                 (tmpmol->get_frame(0)->pos)[ipdb*3+1],
                                 (tmpmol->get_frame(0)->pos)[ipdb*3+2]);
      ipos++;
      if (!use_pdb_field && current_index == indices.end())
        break;
    }

    if ((ipos < pos.size()) || (current_index != indices.end())) {
      cvm::error("Error: the number of records in the PDB file \""+
                  std::string(pdb_filename)+
                  "\" does not appear to match either the total number of atoms,"+
                  " or the number of coordinates requested at this point("+
                  cvm::to_str(pos.size())+").\n", INPUT_ERROR);
      return COLVARS_ERROR;
    }

  } else {

    // when the PDB contains exactly the number of atoms of the array,
    // ignore the fields and just read coordinates
    for (size_t ia = 0; ia < pos.size(); ia++) {
      pos[ia] = cvm::atom_pos((tmpmol->get_frame(0)->pos)[ia*3],
                               (tmpmol->get_frame(0)->pos)[ia*3+1],
                               (tmpmol->get_frame(0)->pos)[ia*3+2]);
    }
  }

  vmd->molecule_delete(tmpmolid);
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}



int colvarproxy_vmd::load_atoms(char const *pdb_filename,
                                  std::vector<cvm::atom> &atoms,
                                  std::string const &pdb_field_str,
                                  double const pdb_field_value)
{
  if (pdb_field_str.size() == 0) {
    cvm::log("Error: must define which PDB field to use "
                      "in order to define atoms from a PDB file.\n");
    cvm::set_error_bits(INPUT_ERROR);
    return COLVARS_ERROR;
  }

  FileSpec *tmpspec = new FileSpec();
  tmpspec->autobonds = 0;
  int tmpmolid = vmd->molecule_load(-1, pdb_filename, "pdb", tmpspec);
  DrawMolecule *tmpmol = vmd->moleculeList->mol_from_id(tmpmolid);
  delete tmpspec;
  vmd->molecule_make_top(vmdmolid);
  size_t const pdb_natoms = tmpmol->nAtoms;

  e_pdb_field pdb_field_index = pdb_field_str2enum(pdb_field_str);

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
      atom_pdb_field_value = (tmpmol->get_frame(0)->pos)[ipdb*3];
      break;
    case e_pdb_y:
      atom_pdb_field_value = (tmpmol->get_frame(0)->pos)[ipdb*3+1];
      break;
    case e_pdb_z:
      atom_pdb_field_value = (tmpmol->get_frame(0)->pos)[ipdb*3+2];
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

    atoms.push_back(cvm::atom(ipdb+1));
  }

  vmd->molecule_delete(tmpmolid);
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


// atom member functions, VMD specific implementations

cvm::atom::atom(int const &atom_number)
{
  // VMD internal numbering starts from zero
  int const aid(atom_number-1);

  DrawMolecule *vmdmol = ((colvarproxy_vmd *) cvm::proxy)->vmdmol;
  float *masses = vmdmol->mass();

  if (cvm::debug())
    cvm::log("Adding atom "+cvm::to_str(aid+1)+
              " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= vmdmol->nAtoms) ) {
    cvm::error("Error: invalid atom number specified, "+
                      cvm::to_str(atom_number)+"\n");
    return;
  }

  this->id = aid;
  this->mass = masses[aid];
  this->reset_data();
}


// In case of PSF structure, this function's argument "resid" is the non-unique identifier
// TODO: check that the default segment_id of non-PSF topologies is MAIN
cvm::atom::atom(cvm::residue_id const &resid,
                 std::string const     &atom_name,
                 std::string const     &segment_name)
{
  DrawMolecule *vmdmol = ((colvarproxy_vmd *) cvm::proxy)->vmdmol;
  float *masses = vmdmol->mass();

  int aid = -1;
  for (int ir = 0; ir < vmdmol->nResidues; ir++) {
    Residue *vmdres = vmdmol->residue(ir);
    if (vmdres->resid == resid) {
      for (int ia = 0; ia < vmdres->atoms.num(); ia++) {
        int const resaid = vmdres->atoms[ia];
        std::string const sel_segname((vmdmol->segNames).name(vmdmol->atom(resaid)->segnameindex));
        std::string const sel_atom_name((vmdmol->atomNames).name(vmdmol->atom(resaid)->nameindex));
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
    cvm::log("Adding atom \""+
              atom_name+"\" in residue "+
              cvm::to_str(resid)+
              " (index "+cvm::to_str(aid)+
              ") for collective variables calculation.\n");

  if (aid < 0) {
    cvm::error("Error: could not find atom \""+
                      atom_name+"\" in residue "+
                      cvm::to_str(resid)+
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
cvm::atom::atom(cvm::atom const &a)
  : index(a.index), id(a.id), mass(a.mass)
{}


cvm::atom::~atom()
{}

void cvm::atom::read_position()
{
  // read the position directly from the current timestep's memory
  // Note: no prior update should be required (unlike NAMD with GlobalMaster)
  DrawMolecule *vmdmol = ((colvarproxy_vmd *) cvm::proxy)->vmdmol;
  int frame = ((colvarproxy_vmd *) cvm::proxy)->vmdmol_frame;
  float *vmdpos = (vmdmol->get_frame(frame))->pos;
  this->pos = cvm::atom_pos(vmdpos[this->id*3+0],
                             vmdpos[this->id*3+1],
                             vmdpos[this->id*3+2]);
}

void cvm::atom::read_velocity()
{
  // Unavailable, but do not display an error to avoid flooding the output
  return;
}


void cvm::atom::read_system_force()
{
  // Unavailable, but do not display an error to avoid flooding the output
  return;
}


void cvm::atom::apply_force(cvm::rvector const &new_force)
{
  // Unavailable, but do not display an error to avoid flooding the output
  return;
}

