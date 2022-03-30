// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

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


namespace {
  // Keep pointers to relevant runtime objects
  VMDApp *colvars_vmd_ptr = NULL;
  colvarproxy_vmd *colvarproxy_vmd_ptr = NULL;
}


// Copy of declaration from colvarscript.cpp (this alone doesn't warrant a
// new header file)
// COLVARS_TCL is defined when relevant in colvarproxy_tcl.h, included via colvarproxy.h
#ifdef COLVARS_TCL
extern "C" int tcl_run_colvarscript_command(ClientData clientData,
                                            Tcl_Interp *interp_in,
                                            int objc, Tcl_Obj *const objv[]);
#endif


int tcl_colvars(ClientData clientData, Tcl_Interp *interp,
                int objc, Tcl_Obj *const objv[])
{
  // Get pointer to VMD object
  if (colvars_vmd_ptr == NULL) {
    colvars_vmd_ptr = (VMDApp *) clientData;
  }
  return tcl_run_colvarscript_command(clientData, interp, objc, objv);
}


int tcl_colvars_vmd_init(Tcl_Interp *interp, int molid_input)
{
  VMDApp *const vmd = colvars_vmd_ptr;
  int molid = molid_input == -(1<<16) ? vmd->molecule_top() : molid_input;
  if (vmd->molecule_valid_id(molid)) {
    colvarproxy_vmd_ptr = new colvarproxy_vmd(interp, vmd, molid);
    return (cvm::get_error() == COLVARS_OK) ? TCL_OK : TCL_ERROR;
  } else {
    Tcl_SetResult(interp, (char *) "Error: molecule not found.",
                  TCL_STATIC);
    return TCL_ERROR;
  }
}


colvarproxy_vmd::colvarproxy_vmd(Tcl_Interp *interp, VMDApp *v, int molid)
  : vmd(v),
    vmdmolid(molid),
#if defined(VMDTKCON)
    msgColvars("colvars: ", VMDCON_INFO)
#else
    msgColvars("colvars: ")
#endif
{
  version_int = get_version_from_string(COLVARPROXY_VERSION);
  b_simulation_running = false;

  // both fields are taken from data structures already available
  updated_masses_ = updated_charges_ = true;

  // Do we have scripts?
  // For now colvars depend on Tcl, but this may not always be the case
  // in the future
#if defined(VMDTCL)
  have_scripts = true;
  // Need to set this before constructing colvarmodule, which creates colvarscript object
  set_tcl_interp(interp);
#else
  have_scripts = false;
#endif

  colvars = new colvarmodule(this);
  cvm::log("Using VMD interface, version "+
           cvm::to_str(COLVARPROXY_VERSION)+".\n");

  colvars->cite_feature("VMD engine");
  colvars->cite_feature("Colvars-VMD interface (command line)");

  colvars->cv_traj_freq = 0; // I/O will be handled explicitly
  colvars->restart_out_freq = 0;
  cvm::rotation::monitor_crossings = false; // Avoid unnecessary error messages

  // Default to VMD's native unit system, but do not set the units string
  // to preserve the native workflow of VMD / NAMD / LAMMPS-real
  angstrom_value = 1.;
  kcal_mol_value = 1.;

  colvars->setup_input();
  colvars->setup_output();

  // set the same seed as in Measure.C
  vmd_srandom(38572111);

  colvarproxy_vmd::setup();
}


colvarproxy_vmd::~colvarproxy_vmd()
{}


int colvarproxy_vmd::request_deletion()
{
  b_delete_requested = true;
  return COLVARS_OK;
}


int colvarproxy_vmd::setup()
{
  vmdmol = vmd->moleculeList->mol_from_id(vmdmolid);
  if (vmdmol) {
    set_frame(vmdmol->frame());
  } else {
    error("Error: requested molecule ("+cvm::to_str(vmdmolid)+") does not exist.\n");
    return COLVARS_ERROR;
  }

  if (colvars) {
    return colvars->setup();
  }

  return COLVARS_OK;
}


cvm::real colvarproxy_vmd::backend_angstrom_value()
{
  return 1.0;
}


cvm::real colvarproxy_vmd::boltzmann()
{
  return 0.001987191;
}


cvm::real colvarproxy_vmd::temperature()
{
  // TODO define, document and implement a user method to set the value of this
  return 300.0;
}


cvm::real colvarproxy_vmd::dt()
{
  // TODO define, document and implement a user method to set the value of this
  return 1.0;
}


cvm::real colvarproxy_vmd::rand_gaussian()
{
  return vmd_random_gaussian();
}


int colvarproxy_vmd::set_unit_system(std::string const &units_in, bool check_only)
{
  // if check_only is specified, just test for compatibility
  // cvolvarmodule does that if new units are requested while colvars are already defined
  if (check_only) {
    if ((units != "" && units_in != units) || (units == "" && units_in != "real")) {
      cvm::error("Specified unit system \"" + units_in + "\" is incompatible with previous setting \""
                  + units + "\".\nReset the Colvars Module or delete all variables to change the unit.\n");
      return COLVARS_ERROR;
    } else {
      return COLVARS_OK;
    }
  }

  if (units_in == "real") {
    angstrom_value = 1.;
    kcal_mol_value = 1.;
  } else if (units_in == "metal") {
    angstrom_value = 1.;
    kcal_mol_value = 0.0433641017; // eV
    // inverse of LAMMPS value is 1/23.060549 = .043364102
  } else if (units_in == "electron") {
    angstrom_value = 1.88972612;    // Bohr
    kcal_mol_value = 0.00159360144; // Hartree
  } else if (units_in == "gromacs") {
    angstrom_value = 0.1;    // nm
    kcal_mol_value = 4.184;  // kJ/mol
  } else {
    cvm::error("Unknown unit system specified: \"" + units_in + "\". Supported are real, metal, electron, and gromacs.\n");
    return COLVARS_ERROR;
  }

  units = units_in;
  return COLVARS_OK;
}


int colvarproxy_vmd::update_input()
{
  colvarproxy::update_input();

  int error_code = COLVARS_OK;

  // Check that our parent molecule still exists
  if (vmd->moleculeList->mol_from_id(vmdmolid) == NULL) {
    error("Error: requested molecule ("+cvm::to_str(vmdmolid)+") does not exist.\n");
    return COLVARS_ERROR;
  }
  error_code |= update_atomic_properties();

  size_t i;
  // We're not applying any forces but they can be tracked through [cv getatomappliedforces]
  // Clear before updating Module
  for (i = 0; i < atoms_new_colvar_forces.size(); i++) {
    atoms_new_colvar_forces[i].reset();
  }

  // Do we still have a valid frame?
  if (error_code || vmdmol->get_frame(vmdmol_frame) == NULL) {
    error_code |= COLVARS_NO_SUCH_FRAME;
    return error_code;
  }

  // copy positions in the internal arrays
  float *vmdpos = (vmdmol->get_frame(vmdmol_frame))->pos;
  for (i = 0; i < atoms_positions.size(); i++) {
    atoms_positions[i] = cvm::atom_pos(angstrom_to_internal(vmdpos[atoms_ids[i]*3+0]),
                                       angstrom_to_internal(vmdpos[atoms_ids[i]*3+1]),
                                       angstrom_to_internal(vmdpos[atoms_ids[i]*3+2]));
  }


  Timestep const *ts = vmdmol->get_frame(vmdmol_frame);
  {
    // Get lattice vectors
    float A[3];
    float B[3];
    float C[3];
    ts->get_transform_vectors(A, B, C);
    unit_cell_x.set(angstrom_to_internal(A[0]), angstrom_to_internal(A[1]), angstrom_to_internal(A[2]));
    unit_cell_y.set(angstrom_to_internal(B[0]), angstrom_to_internal(B[1]), angstrom_to_internal(B[2]));
    unit_cell_z.set(angstrom_to_internal(C[0]), angstrom_to_internal(C[1]), angstrom_to_internal(C[2]));
  }

  if (ts->a_length + ts->b_length + ts->c_length < 1.0e-6) {
    boundaries_type = boundaries_non_periodic;
    reset_pbc_lattice();
  } else if ((ts->a_length > 1.0e-6) &&
             (ts->b_length > 1.0e-6) &&
             (ts->c_length > 1.0e-6)) {
    if (((ts->alpha-90.0)*(ts->alpha-90.0)) +
        ((ts->beta-90.0)*(ts->beta-90.0)) +
        ((ts->gamma-90.0)*(ts->gamma-90.0)) < 1.0e-6) {
      boundaries_type = boundaries_pbc_ortho;
    } else {
      boundaries_type = boundaries_pbc_triclinic;
    }
    colvarproxy_system::update_pbc_lattice();
  } else {
    boundaries_type = boundaries_unsupported;
  }

  return error_code;
}


void colvarproxy_vmd::add_energy(cvm::real energy)
{
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
  add_error_msg(message);
  log(message);
}


int colvarproxy_vmd::get_molid(int &molid)
{
  molid = vmdmolid;
  return COLVARS_OK;
}


int colvarproxy_vmd::get_frame(long int &f)
{
  f = vmdmol_frame;
  return COLVARS_OK;
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


void colvarproxy_vmd::init_tcl_pointers()
{
#ifdef VMDTCL
  // Do nothing (when constructing colvarproxy), this will be set by colvarproxy_vmd constructor
#else
  colvarproxy::init_tcl_pointers(); // Create dedicated interpreter
#endif
}


// Callback functions

#ifdef VMDTCL

int colvarproxy_vmd::run_force_callback()
{
  return colvarproxy::tcl_run_force_callback();
}

int colvarproxy_vmd::run_colvar_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         colvarvalue &value)
{
  return colvarproxy::tcl_run_colvar_callback(name, cvc_values, value);
}

int colvarproxy_vmd::run_colvar_gradient_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
  return colvarproxy::tcl_run_colvar_gradient_callback(name, cvc_values,
                                                       gradient);
}
#endif


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

      pos[ipos] = cvm::atom_pos(angstrom_to_internal((tmpmol->get_frame(0)->pos)[ipdb*3]),
                                angstrom_to_internal((tmpmol->get_frame(0)->pos)[ipdb*3+1]),
                                angstrom_to_internal((tmpmol->get_frame(0)->pos)[ipdb*3+2]));
      ipos++;
      if (!use_pdb_field && current_index == indices.end())
        break;
    }

    if (ipos < pos.size() || (!use_pdb_field && current_index != indices.end())) {
      size_t n_requested = use_pdb_field ? pos.size() : indices.size();
      cvm::error("Error: number of matching records in the PDB file \""+
                 std::string(pdb_filename)+"\" ("+cvm::to_str(ipos)+
                 ") does not match the number of requested coordinates ("+
                 cvm::to_str(n_requested)+").\n", INPUT_ERROR);
      return COLVARS_ERROR;
    }
  } else {

    // when the PDB contains exactly the number of atoms of the array,
    // ignore the fields and just read coordinates
    for (size_t ia = 0; ia < pos.size(); ia++) {
      pos[ia] = cvm::atom_pos(angstrom_to_internal((tmpmol->get_frame(0)->pos)[ia*3]),
                              angstrom_to_internal((tmpmol->get_frame(0)->pos)[ia*3+1]),
                              angstrom_to_internal((tmpmol->get_frame(0)->pos)[ia*3+2]));
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
  delete tmpspec;

  if (tmpmolid < 0) {
    cvm::error("Error: VMD could not read file \""+std::string(pdb_filename)+"\".\n",
               FILE_ERROR);
    return COLVARS_ERROR;
  }
  DrawMolecule *tmpmol = vmd->moleculeList->mol_from_id(tmpmolid);

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
    cvm::log("Adding atom \""+
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


int colvarproxy_vmd::init_volmap_by_id(int volmap_id)
{
  for (size_t i = 0; i < volmaps_ids.size(); i++) {
    if (volmaps_ids[i] == volmap_id) {
      // this map has already been requested
      volmaps_ncopies[i] += 1;
      return i;
    }
  }

  int index = -1;

  int error_code = check_volmap_by_id(volmap_id);
  if (error_code == COLVARS_OK) {
    index = add_volmap_slot(volmap_id);
  }

  return index;
}


int colvarproxy_vmd::check_volmap_by_id(int volmap_id)
{
  if ((volmap_id < 0) || (volmap_id >= vmdmol->num_volume_data())) {
    return cvm::error("Error: invalid numeric ID ("+cvm::to_str(volmap_id)+
                      ") for map.\n", INPUT_ERROR);
  }
  return COLVARS_OK;
}


void colvarproxy_vmd::clear_volmap(int index)
{
  colvarproxy::clear_volmap(index);
}


template<int flags>
void colvarproxy_vmd::compute_voldata(VolumetricData const *voldata,
                                      cvm::atom_iter atom_begin,
                                      cvm::atom_iter atom_end,
                                      cvm::real *value,
                                      cvm::real *atom_field)
{
  int i = 0;
  float coord[3], voxcoord[3], grad[3];
  cvm::rvector dV(0.0);
  cvm::atom_iter ai = atom_begin;
  cvm::atom_pos const origin(0.0, 0.0, 0.0);

  for ( ; ai != atom_end; ai++, i++) {

    // Wrap around the origin
    cvm::rvector const wrapped_pos = position_distance(origin, ai->pos);
    coord[0] = internal_to_angstrom(wrapped_pos.x);
    coord[1] = internal_to_angstrom(wrapped_pos.y);
    coord[2] = internal_to_angstrom(wrapped_pos.z);

    // Get the coordinates on the grid and check the boundaries
    voldata->voxel_coord_from_cartesian_coord(coord, voxcoord, 0);
    int const gx = static_cast<int>(voxcoord[0]);
    int const gy = static_cast<int>(voxcoord[1]);
    int const gz = static_cast<int>(voxcoord[2]);
    int valid_coord = 1;
    if ((gx < 0) || (gx >= voldata->xsize) ||
        (gy < 0) || (gy >= voldata->ysize) ||
        (gz < 0) || (gz >= voldata->zsize)) {
      valid_coord = 0;
    }

    cvm::real const V = valid_coord ?
      static_cast<cvm::real>(voldata->voxel_value_interpolate(voxcoord[0],
                                                              voxcoord[1],
                                                              voxcoord[2])) :
      0.0;

    if (flags & volmap_flag_gradients) {
      if (valid_coord) {
        voldata->voxel_gradient_interpolate(voxcoord, grad);
        // Correct the sign of the gradient
        dV = cvm::rvector(-1.0*grad[0], -1.0*grad[1], -1.0*grad[2]);
      } else {
        dV = cvm::rvector(0.0);
      }
    }

    if (flags & volmap_flag_use_atom_field) {
      *value += V * atom_field[i];
      if (flags & volmap_flag_gradients) {
        ai->grad += atom_field[i] * dV;
      }
    } else {
      *value += V;
      if (flags & volmap_flag_gradients) {
        ai->grad += dV;
      }
    }
  }
}


int colvarproxy_vmd::compute_volmap(int flags,
                                    int volmap_id,
                                    cvm::atom_iter atom_begin,
                                    cvm::atom_iter atom_end,
                                    cvm::real *value,
                                    cvm::real *atom_field)
{
  int error_code = COLVARS_OK;
  VolumetricData const *voldata = vmdmol->get_volume_data(volmap_id);
  if (voldata != NULL) {

    if (flags & volmap_flag_gradients) {

      if (flags & volmap_flag_use_atom_field) {
        int const new_flags = volmap_flag_gradients |
          volmap_flag_use_atom_field;
        compute_voldata<new_flags>(voldata, atom_begin, atom_end,
                                   value, atom_field);
      } else {
        int const new_flags = volmap_flag_gradients;
        compute_voldata<new_flags>(voldata, atom_begin, atom_end,
                                   value, NULL);
      }

    } else {

      if (flags & volmap_flag_use_atom_field) {
        int const new_flags = volmap_flag_use_atom_field;
        compute_voldata<new_flags>(voldata, atom_begin, atom_end,
                                   value, atom_field);
      } else {
        int const new_flags = volmap_flag_null;
        compute_voldata<new_flags>(voldata, atom_begin, atom_end,
                                   value, NULL);
      }
    }
  } else {
    // Error message
    error_code |=
      const_cast<colvarproxy_vmd *>(this)->check_volmap_by_id(volmap_id);
  }
  return error_code;
}
