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
  static std::string tcl_result;
  int script_retval;

  if (proxy != NULL) {

    if ( argc >= 3 ) {
      if (!strcmp(argv[1], "molid")) {
        Tcl_SetResult(interp, (char *) (std::string("Colvars module already created: type \"cv\" for a list of arguments.").c_str()), TCL_STATIC);
        return TCL_ERROR;
      }
    }

    // Clear non-fatal errors from previous commands
    cvm::clear_error();
    proxy->error_output.clear();
    tcl_result.clear();

    script_retval = proxy->script->run(argc, argv);
    // append the error messages from colvarscript to the error messages caught by the proxy
    tcl_result = proxy->error_output + proxy->script->result;
    Tcl_SetResult(interp, (char *) tcl_result.c_str(), TCL_STATIC);

    if (cvm::get_error_bit(DELETE_COLVARS)) {
      delete proxy;
      proxy = NULL;
      return TCL_OK;
    }

    if (cvm::get_error_bit(FATAL_ERROR)) {
      // Fatal error: clean up cvm object and proxy
      delete proxy;
      proxy = NULL;
      return TCL_ERROR;
    }

    if (script_retval == COLVARSCRIPT_OK && !cvm::get_error()) {
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

  total_force_requested = false;

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


int colvarproxy_vmd::setup()
{
  vmdmol = vmd->moleculeList->mol_from_id(vmdmolid);
  if (vmdmol) {
    set_frame(vmdmol->frame());
  } else {
    error("Error: cannot find the molecule requested("+cvm::to_str(vmdmolid)+").\n");
    return COLVARS_ERROR;
  }

  // set the same seed as in Measure.C
  vmd_srandom(38572111);

  if (colvars) {
    return colvars->setup();
  }

  return COLVARS_OK;
}


int colvarproxy_vmd::update_input()
{
  colvarproxy::update_input();

  int error_code = COLVARS_OK;

  error_code |= update_atomic_properties();

  // copy positions in the internal arrays
  float *vmdpos = (vmdmol->get_frame(vmdmol_frame))->pos;
  for (size_t i = 0; i < atoms_positions.size(); i++) {
    atoms_positions[i] = cvm::atom_pos(vmdpos[atoms_ids[i]*3+0],
                                       vmdpos[atoms_ids[i]*3+1],
                                       vmdpos[atoms_ids[i]*3+2]);
  }

  return error_code;
}


int colvarproxy_vmd::update_atomic_properties()
{
  float const *masses = vmdmol->mass();
  float const *charges = vmdmol->charge();

  int error_code = COLVARS_OK;

  if (masses == NULL) {
    error("Error: masses are undefined for the molecule being used.\n");
    error_code |= BUG_ERROR;
  } else {
    for (size_t i = 0; i < atoms_ids.size(); i++) {
      atoms_masses[i]  = masses[atoms_ids[i]];
    }
  }

  if (charges == NULL) {
    error("Error: charges are undefined for the molecule being used.\n");
    error_code |= BUG_ERROR;
  } else {
    for (size_t i = 0; i < atoms_ids.size(); i++) {
      atoms_charges[i] = charges[atoms_ids[i]];
    }
  }

  return error_code;
}


void colvarproxy_vmd::request_total_force(bool yesno)
{
  if ((yesno == true) && (total_force_requested == false)) {
    cvm::log("Warning: a bias requested total forces, which are undefined in VMD.  "
             "This is only meaningful when analyzing a simulation where these were used, "
             "provided that a state file is loaded.\n");
  }
  total_force_requested = yesno;
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
  error_output += message;
  log(message);
  if (!cvm::debug()) {
    log("If this error message is unclear, "
        "try recompiling VMD with -DCOLVARS_DEBUG.\n");
  }
}


void colvarproxy_vmd::fatal_error(std::string const &message)
{
  // Fatal error bit is already set, will be handled
  // by tcl_colvars() before handing control back to VMD
  error(message);
}


void colvarproxy_vmd::exit(std::string const &message)
{
  error("Error: requested VMD shutdown.\n");
}


int colvarproxy_vmd::set_frame(long int f)
{
  if (vmdmol->get_frame(f) != NULL) {
    vmdmol_frame = f;
    colvars->it = f;
    update_input();
    return COLVARS_OK;
  } else {
    return COLVARS_NO_SUCH_FRAME;
  }
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

  e_pdb_field pdb_field_index = e_pdb_none;
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
                                cvm::atom_group &atoms,
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

    atoms.add_atom(cvm::atom(ipdb+1));
  }

  vmd->molecule_delete(tmpmolid);
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}



int colvarproxy_vmd::check_atom_id(int atom_number)
{
  // VMD internal numbering starts from zero
  int const aid(atom_number-1);

  if (cvm::debug())
    cvm::log("Adding atom "+cvm::to_str(aid+1)+
             " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= vmdmol->nAtoms) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n");
    return INPUT_ERROR;
  }

  return aid;
}


int colvarproxy_vmd::init_atom(int atom_number)
{
  // save time by checking first whether this atom has been requested before
  // (this is more common than a non-valid atom number)
  int aid = (atom_number-1);

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_ncopies[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);

  if (aid < 0) {
    return INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);

  float const *masses = vmdmol->mass();
  float const *charges = vmdmol->charge();
  atoms_masses[index] = masses[aid];
  atoms_charges[index] = charges[aid];

  return index;
}


int colvarproxy_vmd::check_atom_id(cvm::residue_id const &resid,
                                   std::string const     &atom_name,
                                   std::string const     &segment_id)
{
  int aid = -1;
  for (int ir = 0; ir < vmdmol->nResidues; ir++) {
    Residue *vmdres = vmdmol->residue(ir);
    if (vmdres->resid == resid) {
      for (int ia = 0; ia < vmdres->atoms.num(); ia++) {
        int const resaid = vmdres->atoms[ia];
        std::string const sel_segname((vmdmol->segNames).name(vmdmol->atom(resaid)->segnameindex));
        std::string const sel_atom_name((vmdmol->atomNames).name(vmdmol->atom(resaid)->nameindex));
        if ( ((segment_id.size() == 0) || (segment_id == sel_segname)) &&
             (atom_name == sel_atom_name) ) {
          aid = resaid;
          break;
        }
      }
    }
    if (aid >= 0) break;
  }


  if (aid < 0) {
    cvm::error("Error: could not find atom \""+
               atom_name+"\" in residue "+
               cvm::to_str(resid)+
               ( (segment_id.size()) ?
                 (", segment \""+segment_id+"\"") :
                 ("") )+
               "\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  return aid;
}


int colvarproxy_vmd::init_atom(cvm::residue_id const &resid,
                               std::string const     &atom_name,
                               std::string const     &segment_id)
{
  int const aid = check_atom_id(resid, atom_name, segment_id);

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_ncopies[i] += 1;
      return i;
    }
  }

  if (cvm::debug())
    log("Adding atom \""+
        atom_name+"\" in residue "+
        cvm::to_str(resid)+
        " (index "+cvm::to_str(aid)+
        ") for collective variables calculation.\n");

  int const index = add_atom_slot(aid);

  float const *masses = vmdmol->mass();
  float const *charges = vmdmol->charge();
  atoms_masses[index] = masses[aid];
  atoms_charges[index] = charges[aid];

  return index;
}

