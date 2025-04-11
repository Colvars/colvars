// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <errno.h>

#include "common.h"
#include "fstream_namd.h"
#include "Debug.h"
#include "BackEnd.h"
#include "InfoStream.h"
#include "Node.h"
#include "Molecule.h"
#include "GridForceGrid.h"
#include "GridForceGrid.inl"
#include "PDB.h"
#include "PDBData.h"
#include "ReductionMgr.h"
#include "ScriptTcl.h"
#include "NamdState.h"
#include "Controller.h"
#include "PatchData.h"

#ifdef NAMD_TCL
#include <tcl.h>
#endif

// For replica exchange
#include "converse.h"
#include "DataExchanger.h"

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "colvarproxy_namd.h"
#include "colvarscript.h"


colvarproxy_namd::colvarproxy_namd()
{
  engine_name_ = "NAMD";
#if CMK_SMP && USE_CKLOOP
  charm_lock_state = CmiCreateLock();
#endif

  version_int = get_version_from_string(COLVARPROXY_VERSION);
#if CMK_TRACE_ENABLED
  if ( 0 == CkMyPe() ) {
    traceRegisterUserEvent("GM COLVAR item", GLOBAL_MASTER_CKLOOP_CALC_ITEM);
    traceRegisterUserEvent("GM COLVAR bias", GLOBAL_MASTER_CKLOOP_CALC_BIASES );
    traceRegisterUserEvent("GM COLVAR scripted bias", GLOBAL_MASTER_CKLOOP_CALC_SCRIPTED_BIASES );
  }
#endif
  first_timestep = true;
  requestTotalForce(total_force_requested);

  boltzmann_ = 0.001987191;

  angstrom_value_ = 1.;

  // initialize pointers to NAMD configuration data
  simparams = Node::Object()->simParameters;

  if (cvm::debug())
    iout << "Info: initializing the colvars proxy object.\n" << endi;

  // find the configuration file, if provided
  StringList *config = Node::Object()->configList->find("colvarsConfig");

  // find the input state file
  StringList *input_restart = Node::Object()->configList->find("colvarsInput");
  colvarproxy_io::set_input_prefix(input_restart ? input_restart->data : "");

  update_target_temperature();
  set_integration_timestep(simparams->dt);

  random = Random(simparams->randomSeed);

  // both fields are taken from data structures already available
  updated_masses_ = updated_charges_ = true;

  // Take the output prefixes from the NAMD input
  colvarproxy_io::set_output_prefix(std::string(simparams->outputFilename));
  colvarproxy_io::set_restart_output_prefix(std::string(simparams->restartFilename));
  colvarproxy_io::set_default_restart_frequency(simparams->restartFrequency);

  if (simparams->accelMDOn) {
    accelMDOn = true;
  } else {
    accelMDOn = false;
  }
  amd_weight_factor = 1.0;

  // check if it is possible to save output configuration
  if ((!output_prefix_str.size()) && (!restart_output_prefix_str.size())) {
    error("Error: neither the final output state file or "
          "the output restart file could be defined, exiting.\n");
  }

  init_atoms_map();

  // initialize module: this object will be the communication proxy
  colvars = new colvarmodule(this);

  cvm::log("Using NAMD interface, version "+
           cvm::to_str(COLVARPROXY_VERSION)+".\n");
  colvars->cite_feature("NAMD engine");
  colvars->cite_feature("Colvars-NAMD interface");

  errno = 0;
  for ( ; config; config = config->next ) {
    add_config("configfile", config->data);
  }

  // Trigger immediate initialization of the module
  colvarproxy::parse_module_config();
  colvarproxy_namd::setup();
  colvars->update_engine_parameters();
  colvars->setup_input();
  colvars->setup_output();

  // save to Node for Tcl script access
  Node::Object()->colvars = colvars;

#ifdef NAMD_TCL
  have_scripts = true;
#else
  have_scripts = false;
#endif

  if (simparams->firstTimestep != 0) {
    colvars->set_initial_step(static_cast<cvm::step_number>(simparams->firstTimestep));
  }

#if !defined (NAMD_UNIFIED_REDUCTION)
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
#endif

  #if defined(NODEGROUP_FORCE_REGISTER) && !defined(NAMD_UNIFIED_REDUCTION)
  CProxy_PatchData cpdata(CkpvAccess(BOCclass_group).patchData);
  PatchData *patchData = cpdata.ckLocalBranch();
  nodeReduction = patchData->reduction;
  #endif

  if (cvm::debug())
    iout << "Info: done initializing the colvars proxy object.\n" << endi;
}


colvarproxy_namd::~colvarproxy_namd()
{
#if CMK_SMP && USE_CKLOOP
  CmiDestroyLock(charm_lock_state);
#endif
#if !defined (NAMD_UNIFIED_REDUCTION)
  delete reduction;
#endif
}


int colvarproxy_namd::update_target_temperature()
{
  int error_code = COLVARS_OK;
  if (simparams->rescaleFreq > 0) {
    error_code |= set_target_temperature(simparams->rescaleTemp);
  } else if (simparams->reassignFreq > 0) {
    error_code |= set_target_temperature(simparams->reassignTemp);
  } else if (simparams->langevinOn) {
    error_code |= set_target_temperature(simparams->langevinTemp);
  } else if (simparams->tCoupleOn) {
    error_code |= set_target_temperature(simparams->tCoupleTemp);
  } else if (simparams->loweAndersenOn) {
    error_code |= set_target_temperature(simparams->loweAndersenTemp);
  } else if (simparams->stochRescaleOn) {
    error_code |= set_target_temperature(simparams->stochRescaleTemp);
  } else {
    error_code |= set_target_temperature(0.0);
  }
  return error_code;
}



void colvarproxy_namd::init_atoms_map()
{
  size_t const n_all_atoms = Node::Object()->molecule->numAtoms;
  atoms_map.assign(n_all_atoms, -1);
}


int colvarproxy_namd::update_atoms_map(AtomIDList::const_iterator begin,
                                       AtomIDList::const_iterator end)
{
  if (atoms_map.size() != Node::Object()->molecule->numAtoms) {
    init_atoms_map();
  }

  for (AtomIDList::const_iterator a_i = begin; a_i != end; a_i++) {

    if (atoms_map[*a_i] >= 0) continue;

    for (size_t i = 0; i < atoms_ids.size(); i++) {
      if (atoms_ids[i] == *a_i) {
        atoms_map[*a_i] = i;
        break;
      }
    }

    if (atoms_map[*a_i] < 0) {
      // this atom is probably managed by another GlobalMaster:
      // add it here anyway to avoid having to test for array boundaries at each step
      int const index = add_atom_slot(*a_i);
      atoms_map[*a_i] = index;
      modifyRequestedAtoms().add(*a_i);
      update_atom_properties(index);
    }
  }

  if (cvm::debug()) {
    cvm::log("atoms_map = "+cvm::to_str(atoms_map)+".\n");
  }

  return COLVARS_OK;
}


int colvarproxy_namd::setup()
{
  int error_code = colvarproxy::setup();

  if (colvars->size() == 0) {
    // Module is empty, nothing to do
    return COLVARS_OK;
  }

  log("Updating NAMD interface:\n");

  errno = 0;

  if (simparams->wrapAll) {
    log("Warning: enabling wrapAll can lead to inconsistent results "
        "for Colvars calculations: please disable wrapAll, "
        "as is the default option in NAMD.\n");
  }

  log("updating atomic data ("+cvm::to_str(atoms_ids.size())+" atoms).\n");

  size_t i;
  for (i = 0; i < atoms_ids.size(); i++) {
    update_atom_properties(i);

    // zero out mutable arrays
    atoms_positions[i] = cvm::rvector(0.0, 0.0, 0.0);
    atoms_total_forces[i] = cvm::rvector(0.0, 0.0, 0.0);
    atoms_new_colvar_forces[i] = cvm::rvector(0.0, 0.0, 0.0);
  }

  size_t n_group_atoms = 0;
  for (int ig = 0; ig < modifyRequestedGroups().size(); ig++) {
    n_group_atoms += modifyRequestedGroups()[ig].size();
  }

  log("updating group data ("+cvm::to_str(atom_groups_ids.size())+
      " scalable groups, "+
      cvm::to_str(n_group_atoms)+" atoms in total).\n");

  // Note: groupMassBegin, groupMassEnd may be used here, but they won't work for charges
  for (int ig = 0; ig < modifyRequestedGroups().size(); ig++) {

    // update mass and charge
    update_group_properties(ig);

    atom_groups_coms[ig] = cvm::rvector(0.0, 0.0, 0.0);
    atom_groups_total_forces[ig] = cvm::rvector(0.0, 0.0, 0.0);
    atom_groups_new_colvar_forces[ig] = cvm::rvector(0.0, 0.0, 0.0);
  }

#if NAMD_VERSION_NUMBER >= 34471681
  log("updating grid object data ("+cvm::to_str(volmaps_ids.size())+
      " grid objects in total).\n");
  for (int imap = 0; imap < modifyGridObjForces().size(); imap++) {
    volmaps_new_colvar_forces[imap] = 0.0;
  }
#endif

  size_t const new_features_hash =
    std::hash<std::string>{}(colvars->feature_report(0));
  if (new_features_hash != features_hash) {
    // Nag only once, there may be many run commands
    log(std::string("\n")+colvars->feature_report(0)+std::string("\n"));
    features_hash = new_features_hash;
  }

  update_target_temperature();
  log("updating target temperature (T = "+
      cvm::to_str(target_temperature())+" K).\n");

  // Note: not needed currently, but may be in the future if NAMD allows
  // redefining the timestep
  set_integration_timestep(simparams->dt);

  return error_code;
}


int colvarproxy_namd::reset()
{
  if (cvm::debug()) {
    cvm::log("colvarproxy_namd::reset()\n");
  }

  int error_code = COLVARS_OK;

  // Unrequest all positions, total forces, etc from NAMD
  modifyRequestedAtoms().clear();
  modifyForcedAtoms().clear();
  modifyAppliedForces().clear();

  modifyRequestedGroups().clear();
  modifyGroupForces().clear();

#if NAMD_VERSION_NUMBER >= 34471681
  modifyRequestedGridObjects().clear();
  modifyGridObjForces().clear();
#endif

  requestTotalForce(false);

  atoms_map.clear();

  // Clear internal atomic data
  error_code |= colvarproxy::reset();

  return error_code;
}


void colvarproxy_namd::calculate()
{
  errno = 0;

  if (first_timestep) {

    // First run after the proxy is constructed

    colvarproxy_namd::setup();
    colvars->update_engine_parameters();
    colvars->setup_input();
    colvars->setup_output();
    // Controller is only available after full startup phase, so now
    controller = &Node::Object()->state->getController();

    first_timestep = false;

  } else {

    // Use the time step number inherited from GlobalMaster
    if ( step - previous_NAMD_step == 1 ) {
      colvars->it++;
      b_simulation_continuing = false;
    } else {

      // Cases covered by this condition:
      // - run 0
      // - beginning of a new run statement
      // The internal counter is not incremented, and the objects are made
      // aware of this via the following flag
      b_simulation_continuing = true;

      // Update NAMD output and restart prefixes
      colvarproxy_io::set_output_prefix(std::string(simparams->outputFilename));
      colvarproxy_io::set_restart_output_prefix(std::string(simparams->restartFilename));
      colvarproxy_io::set_default_restart_frequency(simparams->restartFrequency);
      colvars->setup_output();

    }
  }

  previous_NAMD_step = step;
  if (accelMDOn) update_accelMD_info();

  {
    Vector const a = lattice->a();
    Vector const b = lattice->b();
    Vector const c = lattice->c();
    unit_cell_x.set(a.x, a.y, a.z);
    unit_cell_y.set(b.x, b.y, c.z);
    unit_cell_z.set(c.x, c.y, c.z);
  }

  if (!lattice->a_p() && !lattice->b_p() && !lattice->c_p()) {
    boundaries_type = boundaries_non_periodic;
    reset_pbc_lattice();
  } else if (lattice->a_p() && lattice->b_p() && lattice->c_p()) {
    if (lattice->orthogonal()) {
      boundaries_type = boundaries_pbc_ortho;
    } else {
      boundaries_type = boundaries_pbc_triclinic;
    }
    colvarproxy_system::update_pbc_lattice();
  } else {
    boundaries_type = boundaries_unsupported;
  }

  if (cvm::debug()) {
    cvm::log(std::string(cvm::line_marker)+
             "colvarproxy_namd, step no. "+cvm::to_str(colvars->it)+"\n"+
             "Updating atomic data arrays.\n");
  }

  // must delete the forces applied at the previous step: we can do
  // that because they have already been used and copied to other
  // memory locations
  modifyForcedAtoms().clear();
  modifyAppliedForces().clear();

  // If new atomic positions or forces have been requested by other
  // GlobalMaster objects, add these to the atom map as well
  size_t const n_all_atoms = Node::Object()->molecule->numAtoms;
  if ( (atoms_map.size() != n_all_atoms) ||
       (int(atoms_ids.size()) < (getAtomIdEnd() - getAtomIdBegin())) ||
       (int(atoms_ids.size()) < (getForceIdEnd() - getForceIdBegin())) ) {
    update_atoms_map(getAtomIdBegin(), getAtomIdEnd());
    update_atoms_map(getForceIdBegin(), getForceIdEnd());
  }

  // prepare local arrays
  for (size_t i = 0; i < atoms_ids.size(); i++) {
    atoms_positions[i] = cvm::rvector(0.0, 0.0, 0.0);
    atoms_total_forces[i] = cvm::rvector(0.0, 0.0, 0.0);
    atoms_new_colvar_forces[i] = cvm::rvector(0.0, 0.0, 0.0);
  }

  for (size_t i = 0; i < atom_groups_ids.size(); i++) {
    atom_groups_total_forces[i] = cvm::rvector(0.0, 0.0, 0.0);
    atom_groups_new_colvar_forces[i] = cvm::rvector(0.0, 0.0, 0.0);
  }

#if NAMD_VERSION_NUMBER >= 34471681
  for (int imap = 0; imap < volmaps_ids.size(); imap++) {
    volmaps_new_colvar_forces[imap] = 0.0;
  }
#endif

  {
    if (cvm::debug()) {
      cvm::log("Updating positions arrays.\n");
    }
    size_t n_positions = 0;
    AtomIDList::const_iterator a_i = getAtomIdBegin();
    AtomIDList::const_iterator a_e = getAtomIdEnd();
    PositionList::const_iterator p_i = getAtomPositionBegin();

    for ( ; a_i != a_e; ++a_i, ++p_i ) {
      atoms_positions[atoms_map[*a_i]] = cvm::rvector((*p_i).x, (*p_i).y, (*p_i).z);
      n_positions++;
    }

    // The following had to be relaxed because some atoms may be forced without their position being requested
    // if (n_positions < atoms_ids.size()) {
    //   cvm::error("Error: did not receive the positions of all atoms.\n", COLVARS_BUG_ERROR);
    // }
  }

  if (total_force_requested && cvm::step_relative() > 0) {

    // sort the force arrays the previous step
    // (but only do so if there *is* a previous step!)

    {
      if (cvm::debug()) {
        cvm::log("Updating total forces arrays.\n");
      }
      size_t n_total_forces = 0;
      AtomIDList::const_iterator a_i = getForceIdBegin();
      AtomIDList::const_iterator a_e = getForceIdEnd();
      ForceList::const_iterator f_i = getTotalForce();

      for ( ; a_i != a_e; ++a_i, ++f_i ) {
        atoms_total_forces[atoms_map[*a_i]] = cvm::rvector((*f_i).x, (*f_i).y, (*f_i).z);
        n_total_forces++;
      }

      if ( (! b_simulation_continuing) &&
           (n_total_forces < atoms_ids.size()) ) {
        cvm::error("Error: total forces were requested, but total forces "
                   "were not received for all atoms.\n"
                   "The most probable cause is combination of energy "
                   "minimization with a biasing method that requires MD (e.g. ABF).\n"
                   "Always run minimization and ABF separately.", COLVARS_INPUT_ERROR);
      }
    }

    {
      if (cvm::debug()) {
        cvm::log("Updating group total forces arrays.\n");
      }
      ForceList::const_iterator f_i = getGroupTotalForceBegin();
      ForceList::const_iterator f_e = getGroupTotalForceEnd();
      size_t i = 0;
      if ( (! b_simulation_continuing) &&
           ((f_e - f_i) != ((int) atom_groups_ids.size())) ) {
        cvm::error("Error: total forces were requested for scalable groups, "
                   "but they are not in the same number from the number of groups.\n"
                   "The most probable cause is combination of energy "
                   "minimization with a biasing method that requires MD (e.g. ABF).\n"
                   "Always run minimization and ABF separately.", COLVARS_INPUT_ERROR);
      }
      for ( ; f_i != f_e; f_i++, i++) {
        atom_groups_total_forces[i] = cvm::rvector((*f_i).x, (*f_i).y, (*f_i).z);
      }
    }
  }

  {
    if (cvm::debug()) {
      cvm::log("Updating group positions arrays.\n");
    }
    // update group data (only coms available so far)
    size_t ig;
    // note: getGroupMassBegin() could be used here, but masses and charges
    // have already been calculated from the last call to setup()
    PositionList::const_iterator gp_i = getGroupPositionBegin();
    for (ig = 0; gp_i != getGroupPositionEnd(); gp_i++, ig++) {
      atom_groups_coms[ig] = cvm::rvector(gp_i->x, gp_i->y, gp_i->z);
    }
  }

#if NAMD_VERSION_NUMBER >= 34471681
  {
    if (cvm::debug()) {
      log("Updating grid objects.\n");
    }
    // Using a simple nested loop: there probably won't be so many maps that
    // this becomes performance-limiting
    IntList::const_iterator goi_i = getGridObjIndexBegin();
    BigRealList::const_iterator gov_i = getGridObjValueBegin();
    for ( ; gov_i != getGridObjValueEnd(); goi_i++, gov_i++) {
      for (size_t imap = 0; imap < volmaps_ids.size(); imap++) {
        if (volmaps_ids[imap] == *goi_i) {
          volmaps_values[imap] = *gov_i;
          break;
        }
      }
    }
  }
#endif

  if (cvm::debug()) {
    print_input_atomic_data();
  }

  // call the collective variable module
  if (colvars->calc() != COLVARS_OK) {
    cvm::error("Error in the collective variables module.\n", COLVARS_ERROR);
  }

  if (cvm::debug()) {
    print_output_atomic_data();
  }

  // communicate all forces to the MD integrator
  for (size_t i = 0; i < atoms_ids.size(); i++) {
    cvm::rvector const &f = atoms_new_colvar_forces[i];
    modifyForcedAtoms().add(atoms_ids[i]);
    modifyAppliedForces().add(Vector(f.x, f.y, f.z));
  }

  if (atom_groups_new_colvar_forces.size() > 0) {
    modifyGroupForces().resize(requestedGroups().size());
    ForceList::iterator gf_i = modifyGroupForces().begin();
    for (int ig = 0; gf_i != modifyGroupForces().end(); gf_i++, ig++) {
      cvm::rvector const &f = atom_groups_new_colvar_forces[ig];
      *gf_i = Vector(f.x, f.y, f.z);
    }
  }

#if NAMD_VERSION_NUMBER >= 34471681
  if (volmaps_new_colvar_forces.size() > 0) {
    modifyGridObjForces().resize(requestedGridObjs().size());
    modifyGridObjForces().setall(0.0);
    IntList::const_iterator goi_i = getGridObjIndexBegin();
    BigRealList::iterator gof_i = modifyGridObjForces().begin();
    for ( ; goi_i != getGridObjIndexEnd(); goi_i++, gof_i++) {
      for (size_t imap = 0; imap < volmaps_ids.size(); imap++) {
        if (volmaps_ids[imap] == *goi_i) {
          *gof_i = volmaps_new_colvar_forces[imap];
          break;
        }
      }
    }
  }
#endif

  // send MISC energy
  #if defined(NODEGROUP_FORCE_REGISTER) && !defined(NAMD_UNIFIED_REDUCTION)
  if(!simparams->CUDASOAintegrate) {
    reduction->submit();
  }
  #else
  #if !defined(NAMD_UNIFIED_REDUCTION)
  reduction->submit();
  #else
  submitReduction();
  #endif
  #endif

  // NAMD does not destruct GlobalMaster objects, so we must remember
  // to write all output files at the end of a run
  if (step == simparams->N) {
    post_run();
  }
}

void colvarproxy_namd::update_accelMD_info() {
  // This aMD factor is from previous step!
  amd_weight_factor = std::exp(controller->accelMDdV / (target_temperature() * boltzmann()));
}


// Callback functions

void colvarproxy_namd::init_tcl_pointers()
{
#ifdef NAMD_TCL
  // Store pointer to NAMD's Tcl interpreter
  set_tcl_interp(Node::Object()->getScript()->interp);
#else
  colvarproxy::init_tcl_pointers(); // Create dedicated interpreter
#endif
}

int colvarproxy_namd::run_force_callback()
{
  return colvarproxy::tcl_run_force_callback();
}

int colvarproxy_namd::run_colvar_callback(
                          std::string const &name,
                          std::vector<const colvarvalue *> const &cvc_values,
                          colvarvalue &value)
{
  return colvarproxy::tcl_run_colvar_callback(name, cvc_values, value);
}

int colvarproxy_namd::run_colvar_gradient_callback(
                          std::string const &name,
                          std::vector<const colvarvalue *> const &cvc_values,
                          std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
  return colvarproxy::tcl_run_colvar_gradient_callback(name, cvc_values,
                                                       gradient);
}


void colvarproxy_namd::add_energy(cvm::real energy)
{
  #if defined(NODEGROUP_FORCE_REGISTER) && !defined(NAMD_UNIFIED_REDUCTION)
  if (simparams->CUDASOAintegrate) {
    nodeReduction->item(REDUCTION_MISC_ENERGY) += energy;
  } else {
    reduction->item(REDUCTION_MISC_ENERGY) += energy;
  }
  #else
  #if !defined(NAMD_UNIFIED_REDUCTION)
  reduction->item(REDUCTION_MISC_ENERGY) += energy;
  #else
  addReductionEnergy(REDUCTION_MISC_ENERGY, energy);
  #endif
  #endif
}

void colvarproxy_namd::request_total_force(bool yesno)
{
  if (cvm::debug()) {
    cvm::log("colvarproxy_namd::request_total_force()\n");
  }
  total_force_requested = yesno;
  requestTotalForce(total_force_requested);
  if (cvm::debug()) {
    cvm::log("colvarproxy_namd::request_total_force() end\n");
  }
}


void colvarproxy_namd::log(std::string const &message)
{
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line))
    iout << "colvars: " << line << "\n";
  iout << endi;
}


void colvarproxy_namd::error(std::string const &message)
{
  log(message);
  switch (cvm::get_error()) {
  case COLVARS_FILE_ERROR:
    errno = EIO; break;
  case COLVARS_NOT_IMPLEMENTED:
    errno = ENOSYS; break;
  case COLVARS_MEMORY_ERROR:
    errno = ENOMEM; break;
  }
  char const *msg = "Error in the collective variables module "
    "(see above for details)";
  if (errno) {
    NAMD_err(msg);
  } else {
    NAMD_die(msg);
  }
}


int colvarproxy_namd::check_atom_id(int atom_number)
{
  // NAMD's internal numbering starts from zero
  int const aid = (atom_number-1);

  if (cvm::debug())
    cvm::log("Adding atom "+cvm::to_str(atom_number)+
        " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= Node::Object()->molecule->numAtoms) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  return aid;
}


int colvarproxy_namd::check_atom_name_selections_available()
{
  return COLVARS_OK;
}


int colvarproxy_namd::init_atom(int atom_number)
{
  // save time by checking first whether this atom has been requested before
  // (this is more common than a non-valid atom number)
  int aid = (atom_number-1);

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_refcount[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);

  if (aid < 0) {
    return COLVARS_INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);
  atoms_map[aid] = index;
  modifyRequestedAtoms().add(aid);
  update_atom_properties(index);
  return index;
}


int colvarproxy_namd::check_atom_id(cvm::residue_id const &residue,
                                    std::string const     &atom_name,
                                    std::string const     &segment_id)
{
  int const aid =
    (segment_id.size() ?
     Node::Object()->molecule->get_atom_from_name(segment_id.c_str(),
                                                  residue,
                                                  atom_name.c_str()) :
     Node::Object()->molecule->get_atom_from_name("MAIN",
                                                  residue,
                                                  atom_name.c_str()));

  if (aid < 0) {
    // get_atom_from_name() has returned an error value
    cvm::error("Error: could not find atom \""+
               atom_name+"\" in residue "+
               cvm::to_str(residue)+
               ( (segment_id != "MAIN") ?
                 (", segment \""+segment_id+"\"") :
                 ("") )+
               "\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  return aid;
}



/// For AMBER topologies, the segment id is automatically set to
/// "MAIN" (the segment id assigned by NAMD's AMBER topology parser),
/// and is therefore optional when an AMBER topology is used
int colvarproxy_namd::init_atom(cvm::residue_id const &residue,
                                std::string const     &atom_name,
                                std::string const     &segment_id)
{
  int const aid = check_atom_id(residue, atom_name, segment_id);

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_refcount[i] += 1;
      return i;
    }
  }

  if (cvm::debug())
    cvm::log("Adding atom \""+
        atom_name+"\" in residue "+
        cvm::to_str(residue)+
        " (index "+cvm::to_str(aid)+
        ") for collective variables calculation.\n");

  int const index = add_atom_slot(aid);
  atoms_map[aid] = index;
  modifyRequestedAtoms().add(aid);
  update_atom_properties(index);
  return index;
}


void colvarproxy_namd::clear_atom(int index)
{
  colvarproxy::clear_atom(index);
  // TODO remove it from GlobalMaster arrays?
}


void colvarproxy_namd::update_atom_properties(int index)
{
  Molecule *mol = Node::Object()->molecule;
  // update mass
  double const mass = mol->atommass(atoms_ids[index]);
  if (mass <= 0.001) {
    this->log("Warning: near-zero mass for atom "+
              cvm::to_str(atoms_ids[index]+1)+
              "; expect unstable dynamics if you apply forces to it.\n");
  }
  atoms_masses[index] = mass;
  // update charge
  atoms_charges[index] = mol->atomcharge(atoms_ids[index]);
}


cvm::rvector colvarproxy_namd::position_distance(cvm::atom_pos const &pos1,
                                                 cvm::atom_pos const &pos2)
  const
{
  Position const p1(pos1.x, pos1.y, pos1.z);
  Position const p2(pos2.x, pos2.y, pos2.z);
  // return p2 - p1
  Vector const d = this->lattice->delta(p2, p1);
  return cvm::rvector(d.x, d.y, d.z);
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
    cvm::error("Error: unsupported PDB field, \""+
               pdb_field_str+"\".\n", COLVARS_INPUT_ERROR);
  }

  return pdb_field;
}


int colvarproxy_namd::load_coords_pdb(char const *pdb_filename,
                                      std::vector<cvm::atom_pos> &pos,
                                      const std::vector<int> &indices,
                                      std::string const &pdb_field_str,
                                      double const pdb_field_value)
{
  if (pdb_field_str.size() == 0 && indices.size() == 0) {
    cvm::error("Bug alert: either PDB field should be defined or list of "
               "atom IDs should be available when loading atom coordinates!\n", COLVARS_BUG_ERROR);
  }

  e_pdb_field pdb_field_index;
  bool const use_pdb_field = (pdb_field_str.size() > 0);
  if (use_pdb_field) {
    pdb_field_index = pdb_field_str2enum(pdb_field_str);
  }

  // next index to be looked up in PDB file (if list is supplied)
  std::vector<int>::const_iterator current_index = indices.begin();

  PDB *pdb = new PDB(pdb_filename);
  size_t const pdb_natoms = pdb->num_atoms();

  if (pos.size() != pdb_natoms) {

    bool const pos_allocated = (pos.size() > 0);

    size_t ipos = 0, ipdb = 0;
    for ( ; ipdb < pdb_natoms; ipdb++) {

      if (use_pdb_field) {
        // PDB field mode: skip atoms with wrong value in PDB field
        double atom_pdb_field_value = 0.0;

        switch (pdb_field_index) {
        case e_pdb_occ:
          atom_pdb_field_value = (pdb->atom(ipdb))->occupancy();
          break;
        case e_pdb_beta:
          atom_pdb_field_value = (pdb->atom(ipdb))->temperaturefactor();
          break;
        case e_pdb_x:
          atom_pdb_field_value = (pdb->atom(ipdb))->xcoor();
          break;
        case e_pdb_y:
          atom_pdb_field_value = (pdb->atom(ipdb))->ycoor();
          break;
        case e_pdb_z:
          atom_pdb_field_value = (pdb->atom(ipdb))->zcoor();
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
        if (((int) ipdb) != *current_index) {
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
                   "more atoms than needed.\n", COLVARS_BUG_ERROR);
      }

      pos[ipos] = cvm::atom_pos((pdb->atom(ipdb))->xcoor(),
                                (pdb->atom(ipdb))->ycoor(),
                                (pdb->atom(ipdb))->zcoor());
      ipos++;
      if (!use_pdb_field && current_index == indices.end())
        break;
    }

    if (ipos < pos.size() || (!use_pdb_field && current_index != indices.end())) {
      size_t n_requested = use_pdb_field ? pos.size() : indices.size();
      cvm::error("Error: number of matching records in the PDB file \""+
                 std::string(pdb_filename)+"\" ("+cvm::to_str(ipos)+
                 ") does not match the number of requested coordinates ("+
                 cvm::to_str(n_requested)+").\n", COLVARS_INPUT_ERROR);
      return COLVARS_ERROR;
    }
  } else {

    // when the PDB contains exactly the number of atoms of the array,
    // ignore the fields and just read coordinates
    for (size_t ia = 0; ia < pos.size(); ia++) {
      pos[ia] = cvm::atom_pos((pdb->atom(ia))->xcoor(),
                              (pdb->atom(ia))->ycoor(),
                              (pdb->atom(ia))->zcoor());
    }
  }

  delete pdb;
  return COLVARS_OK;
}


int colvarproxy_namd::load_atoms_pdb(char const *pdb_filename,
                                     cvm::atom_group &atoms,
                                     std::string const &pdb_field_str,
                                     double const pdb_field_value)
{
  if (pdb_field_str.size() == 0)
    cvm::error("Error: must define which PDB field to use "
               "in order to define atoms from a PDB file.\n", COLVARS_INPUT_ERROR);

  PDB *pdb = new PDB(pdb_filename);
  size_t const pdb_natoms = pdb->num_atoms();

  e_pdb_field pdb_field_index = pdb_field_str2enum(pdb_field_str);

  for (size_t ipdb = 0; ipdb < pdb_natoms; ipdb++) {

    double atom_pdb_field_value = 0.0;

    switch (pdb_field_index) {
    case e_pdb_occ:
      atom_pdb_field_value = (pdb->atom(ipdb))->occupancy();
      break;
    case e_pdb_beta:
      atom_pdb_field_value = (pdb->atom(ipdb))->temperaturefactor();
      break;
    case e_pdb_x:
      atom_pdb_field_value = (pdb->atom(ipdb))->xcoor();
      break;
    case e_pdb_y:
      atom_pdb_field_value = (pdb->atom(ipdb))->ycoor();
      break;
    case e_pdb_z:
      atom_pdb_field_value = (pdb->atom(ipdb))->zcoor();
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

    if (atoms.is_enabled(colvardeps::f_ag_scalable)) {
      atoms.add_atom_id(ipdb);
    } else {
      atoms.add_atom(cvm::atom(ipdb+1));
    }
  }

  delete pdb;
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}

int colvarproxy_namd::load_atoms_pdb(char const *pdb_filename,
                                     cvm::atom_group_soa &atoms,
                                     std::string const &pdb_field_str,
                                     double const pdb_field_value)
{
  if (pdb_field_str.size() == 0)
    cvm::error("Error: must define which PDB field to use "
               "in order to define atoms from a PDB file.\n", COLVARS_INPUT_ERROR);

  PDB *pdb = new PDB(pdb_filename);
  size_t const pdb_natoms = pdb->num_atoms();

  e_pdb_field pdb_field_index = pdb_field_str2enum(pdb_field_str);

  auto modify_atoms = get_atom_modifier();
  for (size_t ipdb = 0; ipdb < pdb_natoms; ipdb++) {

    double atom_pdb_field_value = 0.0;

    switch (pdb_field_index) {
    case e_pdb_occ:
      atom_pdb_field_value = (pdb->atom(ipdb))->occupancy();
      break;
    case e_pdb_beta:
      atom_pdb_field_value = (pdb->atom(ipdb))->temperaturefactor();
      break;
    case e_pdb_x:
      atom_pdb_field_value = (pdb->atom(ipdb))->xcoor();
      break;
    case e_pdb_y:
      atom_pdb_field_value = (pdb->atom(ipdb))->ycoor();
      break;
    case e_pdb_z:
      atom_pdb_field_value = (pdb->atom(ipdb))->zcoor();
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

    if (atoms.is_enabled(colvardeps::f_ag_scalable)) {
      atoms.add_atom_id(ipdb);
    } else {
      // modify_atoms.add_atom(cvm::atom(ipdb+1));
      const int atom_number = ipdb+1;
      const int proxy_index = init_atom(atom_number);
      const int atom_id = get_atom_id(proxy_index);
      const cvm::real atom_mass = get_atom_mass(proxy_index);
      const cvm::real atom_charge = get_atom_charge(proxy_index);
      modify_atoms.add_atom(
        cvm::atom_group_soa::simple_atom{
          .proxy_index = proxy_index,
          .id = atom_id,
          .mass = atom_mass,
          .charge = atom_charge,
          .pos = {0, 0, 0},
          .vel = {0, 0, 0},
          .total_force = {0, 0, 0},
          .grad = {0, 0, 0}});
    }
  }

  delete pdb;
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


std::ostream & colvarproxy_namd::output_stream(std::string const &output_name,
                                               std::string const description)
{
  if (cvm::debug()) {
    cvm::log("Using colvarproxy_namd::output_stream()\n");
  }

  if (!io_available()) {
    cvm::error("Error: trying to access an output file/channel "
               "from the wrong thread.\n", COLVARS_BUG_ERROR);
    return *output_stream_error_;
  }

  if (output_streams_.count(output_name) > 0) {
    return *(output_streams_[output_name]);
  }

  backup_file(output_name.c_str());

  output_streams_[output_name] = new ofstream_namd(output_name.c_str(), std::ios::binary);
  if (! output_streams_[output_name]->good()) {
    cvm::error("Error: cannot write to "+description+" \""+output_name+"\".\n",
               COLVARS_FILE_ERROR);
  }

  return *(output_streams_[output_name]);
}


int colvarproxy_namd::flush_output_stream(std::string const &output_name)
{
  if (!io_available()) {
    return COLVARS_OK;
  }

  if (output_streams_.count(output_name) > 0) {
    (reinterpret_cast<ofstream_namd *>(output_streams_[output_name]))->flush();
    return COLVARS_OK;
  }

  return cvm::error("Error: trying to flush an output file/channel "
                    "that wasn't open.\n", COLVARS_BUG_ERROR);
}


int colvarproxy_namd::flush_output_streams()
{
  if (!io_available()) {
    return COLVARS_OK;
  }

  for (std::map<std::string, std::ostream *>::iterator osi = output_streams_.begin();
       osi != output_streams_.end();
       osi++) {
    (reinterpret_cast<ofstream_namd *>(osi->second))->flush();
  }

  return COLVARS_OK;
}


int colvarproxy_namd::close_output_stream(std::string const &output_name)
{
  if (!io_available()) {
    return cvm::error("Error: trying to access an output file/channel "
                      "from the wrong thread.\n", COLVARS_BUG_ERROR);
  }

  if (output_streams_.count(output_name) > 0) {
    (reinterpret_cast<ofstream_namd *>(output_streams_[output_name]))->close();
    delete output_streams_[output_name];
    output_streams_.erase(output_name);
  }

  return COLVARS_OK;
}


int colvarproxy_namd::close_output_streams()
{
  if (! io_available()) {
    return COLVARS_OK;
  }

  for (std::map<std::string, std::ostream *>::iterator osi = output_streams_.begin();
       osi != output_streams_.end();
       osi++) {
    (reinterpret_cast<ofstream_namd *>(osi->second))->close();
  }
  output_streams_.clear();

  return COLVARS_OK;
}


int colvarproxy_namd::backup_file(char const *filename)
{
  if (std::string(filename).rfind(std::string(".colvars.state")) != std::string::npos) {
    NAMD_backup_file(filename, ".old");
  } else {
    NAMD_backup_file(filename, ".BAK");
  }
  return COLVARS_OK;
}


int colvarproxy_namd::init_atom_group(std::vector<int> const &atoms_ids)
{
  if (cvm::debug())
    cvm::log("Requesting from NAMD a group of size "+cvm::to_str(atoms_ids.size())+
        " for collective variables calculation.\n");

  colvars->cite_feature("Scalable center-of-mass computation (NAMD)");

  // Note: modifyRequestedGroups is supposed to be in sync with the colvarproxy arrays,
  // and to stay that way during a simulation

  // compare this new group to those already allocated inside GlobalMaster
  int ig;
  for (ig = 0; ig < modifyRequestedGroups().size(); ig++) {
    AtomIDList const &namd_group = modifyRequestedGroups()[ig];
    bool b_match = true;

    if (namd_group.size() != ((int) atoms_ids.size())) {
      b_match = false;
    } else {
      int ia;
      for (ia = 0; ia < namd_group.size(); ia++) {
        int const aid = atoms_ids[ia];
        if (namd_group[ia] != aid) {
          b_match = false;
          break;
        }
      }
    }

    if (b_match) {
      if (cvm::debug())
        cvm::log("Group was already added.\n");
      // this group already exists
      atom_groups_refcount[ig] += 1;
      return ig;
    }
  }

  // add this group (note: the argument of add_atom_group_slot() is redundant for NAMD, and provided only for consistency)
  size_t const index = add_atom_group_slot(atom_groups_ids.size());
  modifyRequestedGroups().resize(atom_groups_ids.size());
  // the following is done in calculate()
  // modifyGroupForces().resize(atom_groups_ids.size());
  AtomIDList &namd_group = modifyRequestedGroups()[index];
  namd_group.resize(atoms_ids.size());
  int const n_all_atoms = Node::Object()->molecule->numAtoms;
  for (size_t ia = 0; ia < atoms_ids.size(); ia++) {
    int const aid = atoms_ids[ia];
    if (cvm::debug())
      cvm::log("Adding atom "+cvm::to_str(aid+1)+
          " for collective variables calculation.\n");
    if ( (aid < 0) || (aid >= n_all_atoms) ) {
      cvm::error("Error: invalid atom number specified, "+
                 cvm::to_str(aid+1)+"\n", COLVARS_INPUT_ERROR);
      return -1;
    }
    namd_group[ia] = aid;
  }

  update_group_properties(index);

  if (cvm::debug()) {
    cvm::log("Group has index "+cvm::to_str(index)+"\n");
    cvm::log("modifyRequestedGroups length = "+cvm::to_str(modifyRequestedGroups().size())+
        ", modifyGroupForces length = "+cvm::to_str(modifyGroupForces().size())+"\n");
  }

  return index;
}


void colvarproxy_namd::clear_atom_group(int index)
{
  // do nothing, keep the NAMD arrays in sync with the colvarproxy ones
  colvarproxy::clear_atom_group(index);
}


int colvarproxy_namd::update_group_properties(int index)
{
  AtomIDList const &namd_group = modifyRequestedGroups()[index];
  if (cvm::debug()) {
    cvm::log("Re-calculating total mass and charge for scalable group no. "+cvm::to_str(index+1)+" ("+
             cvm::to_str(namd_group.size())+" atoms).\n");
  }

  cvm::real total_mass = 0.0;
  cvm::real total_charge = 0.0;
  for (int i = 0; i < namd_group.size(); i++) {
    total_mass += Node::Object()->molecule->atommass(namd_group[i]);
    total_charge += Node::Object()->molecule->atomcharge(namd_group[i]);
  }
  atom_groups_masses[index] = total_mass;
  atom_groups_charges[index] = total_charge;

  if (cvm::debug()) {
    cvm::log("total mass = "+cvm::to_str(total_mass)+
             ", total charge = "+cvm::to_str(total_charge)+"\n");
  }

  return COLVARS_OK;
}



int colvarproxy_namd::set_unit_system(std::string const &units_in, bool /*check_only*/)
{
  if (units_in != "real") {
    cvm::error("Error: Specified unit system \"" + units_in + "\" is unsupported in NAMD. Supported units are \"real\" (A, kcal/mol).\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}


#if NAMD_VERSION_NUMBER >= 34471681


int colvarproxy_namd::check_volmaps_available()
{
  return COLVARS_OK;
}


int colvarproxy_namd::init_volmap_by_id(int volmap_id)
{
  for (size_t i = 0; i < volmaps_ids.size(); i++) {
    if (volmaps_ids[i] == volmap_id) {
      // this map has already been requested
      volmaps_refcount[i] += 1;
      return i;
    }
  }

  int error_code = check_volmap_by_id(volmap_id);
  int index = -1;
  if (error_code == COLVARS_OK) {
    index = add_volmap_slot(volmap_id);
    modifyRequestedGridObjects().add(volmap_id);
  }

  return (error_code == COLVARS_OK) ? index : -1;
}


int colvarproxy_namd::init_volmap_by_name(char const *volmap_name)
{
  if (volmap_name == NULL) {
    return cvm::error("Error: no grid object name provided.", COLVARS_INPUT_ERROR);
  }

  int error_code = COLVARS_OK;

  error_code |= check_volmap_by_name(volmap_name);

  int index = -1;
  if (error_code == COLVARS_OK) {

    int volmap_id = simparams->mgridforcelist.index_for_key(volmap_name);

    // Check that the scale factor is correctly set to zero
    Molecule *mol = Node::Object()->molecule;
    GridforceGrid const *grid = mol->get_gridfrc_grid(volmap_id);
    Vector const gfScale = grid->get_scale();
    if ((gfScale.x != 0.0) || (gfScale.y != 0.0) || (gfScale.z != 0.0)) {
      error_code |= cvm::error("Error: GridForce map \""+
                               std::string(volmap_name)+
                               "\" has non-zero scale factors.\n",
                               COLVARS_INPUT_ERROR);
    }

    for (size_t i = 0; i < volmaps_ids.size(); i++) {
      if (volmaps_ids[i] == volmap_id) {
        // this map has already been requested
        volmaps_refcount[i] += 1;
        return i;
      }
    }

    index = add_volmap_slot(volmap_id);
    modifyRequestedGridObjects().add(volmap_id);
  }

  return (error_code == COLVARS_OK) ? index : -1;
}


int colvarproxy_namd::check_volmap_by_id(int volmap_id)
{
  Molecule *mol = Node::Object()->molecule;
  if ((volmap_id < 0) || (volmap_id >= mol->numGridforceGrids)) {
    return cvm::error("Error: invalid numeric ID ("+cvm::to_str(volmap_id)+
                      ") for map.\n", COLVARS_INPUT_ERROR);
  }
  colvars->cite_feature("GridForces volumetric map implementation for NAMD");
  return COLVARS_OK;
}


int colvarproxy_namd::check_volmap_by_name(char const *volmap_name)
{
  if (volmap_name == NULL) {
    return cvm::error("Error: no grid object name provided.", COLVARS_INPUT_ERROR);
  }
  int volmap_id = simparams->mgridforcelist.index_for_key(volmap_name);
  if (volmap_id < 0) {
    return cvm::error("Error: invalid map name \""+std::string(volmap_name)+
                      "\".\n", COLVARS_INPUT_ERROR);
  }
  colvars->cite_feature("GridForces volumetric map implementation for NAMD");
  return COLVARS_OK;
}


void colvarproxy_namd::clear_volmap(int index)
{
  // TODO remove from GlobalMaster
  colvarproxy::clear_volmap(index);
}


int colvarproxy_namd::get_volmap_id_from_name(char const *volmap_name)
{
  int const volmap_id =
    simparams->mgridforcelist.index_for_key(volmap_name);
  if (volmap_id < 0) {
    // Print error
    check_volmap_by_name(volmap_name);
  }
  return volmap_id;
}


template<class T, int flags>
void colvarproxy_namd::GridForceGridLoop(T const *g,
#ifdef COLVARS_USE_SOA
                                         cvm::atom_group_soa* ag,
#else
                                         cvm::atom_iter atom_begin,
                                         cvm::atom_iter atom_end,
#endif // COLVARS_USE_SOA
                                         cvm::real *value,
                                         cvm::real *atom_field)
{
  float V = 0.0f;
  Vector dV(0.0);
#ifdef COLVARS_USE_SOA
  for (size_t i = 0; i < ag->size(); ++i) {
    if (g->compute_VdV(Position(ag->pos_x(i), ag->pos_y(i), ag->pos_z(i)))) {
      // out-of-bounds atom
      V = 0.0f;
      dV = 0.0;
    } else {
      if (flags & volmap_flag_use_atom_field) {
        *value += V * atom_field[i];
        if (flags & volmap_flag_gradients) {
          const cvm::rvector grad = atom_field[i] * cvm::rvector(dV.x, dV.y, dV.z);
          ag->grad_x(i) += grad.x;
          ag->grad_y(i) += grad.y;
          ag->grad_z(i) += grad.z;
        }
      } else {
        *value += V;
        if (flags & volmap_flag_gradients) {
          ag->grad_x(i) += dV.x;
          ag->grad_y(i) += dV.y;
          ag->grad_z(i) += dV.z;
        }
      }
    }
  }
#else
  int i = 0;
  cvm::atom_iter ai = atom_begin;
  for ( ; ai != atom_end; ai++, i++) {
    if (g->compute_VdV(Position(ai->pos.x, ai->pos.y, ai->pos.z), V, dV)) {
      // out-of-bounds atom
      V = 0.0f;
      dV = 0.0;
    } else {
      if (flags & volmap_flag_use_atom_field) {
        *value += V * atom_field[i];
        if (flags & volmap_flag_gradients) {
          ai->grad += atom_field[i] * cvm::rvector(dV.x, dV.y, dV.z);
        }
      } else {
        *value += V;
        if (flags & volmap_flag_gradients) {
          ai->grad += cvm::rvector(dV.x, dV.y, dV.z);
        }
      }
    }
  }
#endif // COLVARS_USE_SOA
}


template<class T>
void colvarproxy_namd::getGridForceGridValue(int flags,
                                             T const *g,
#ifdef COLVARS_USE_SOA
                                             cvm::atom_group_soa* ag,
#else
                                             cvm::atom_iter atom_begin,
                                             cvm::atom_iter atom_end,
#endif // COLVARS_USE_SOA
                                             cvm::real *value,
                                             cvm::real *atom_field)
{
  if (flags & volmap_flag_use_atom_field) {
    int const new_flags = volmap_flag_use_atom_field | volmap_flag_gradients;
#ifdef COLVARS_USE_SOA
    GridForceGridLoop<T, new_flags>(g, ag, value, atom_field);
#else
    GridForceGridLoop<T, new_flags>(g, atom_begin, atom_end,
                                    value, atom_field);
#endif // COLVARS_USE_SOA
  } else {
    int const new_flags = volmap_flag_gradients;
#ifdef COLVARS_USE_SOA
    GridForceGridLoop<T, new_flags>(g, ag, value, atom_field);
#else
    GridForceGridLoop<T, new_flags>(g, atom_begin, atom_end,
                                    value, atom_field);
#endif // COLVARS_USE_SOA
  }
}

#ifdef COLVARS_USE_SOA
int colvarproxy_namd::compute_volmap(int flags,
                                     int volmap_id,
                                     cvm::atom_group_soa* ag,
                                     cvm::real *value,
                                     cvm::real *atom_field)
#else
int colvarproxy_namd::compute_volmap(int flags,
                                     int volmap_id,
                                     cvm::atom_iter atom_begin,
                                     cvm::atom_iter atom_end,
                                     cvm::real *value,
                                     cvm::real *atom_field)
#endif // COLVARS_USE_SOA
{
  Molecule *mol = Node::Object()->molecule;
  GridforceGrid *grid = mol->get_gridfrc_grid(volmap_id);
  // Inheritance is not possible with GridForceGrid's design
  if (grid->get_grid_type() == GridforceGrid::GridforceGridTypeFull) {
    GridforceFullMainGrid *g = dynamic_cast<GridforceFullMainGrid *>(grid);
#ifdef COLVARS_USE_SOA
    getGridForceGridValue<GridforceFullMainGrid>(flags, g, ag,
                                                 value, atom_field);
#else
    getGridForceGridValue<GridforceFullMainGrid>(flags, g, atom_begin, atom_end,
                                                 value, atom_field);
#endif // COLVARS_USE_SOA
  } else if (grid->get_grid_type() == GridforceGrid::GridforceGridTypeLite) {
    GridforceLiteGrid *g = dynamic_cast<GridforceLiteGrid *>(grid);
#ifdef COLVARS_USE_SOA
    getGridForceGridValue<GridforceLiteGrid>(flags, g, ag,
                                             value, atom_field);
#else
    getGridForceGridValue<GridforceLiteGrid>(flags, g, atom_begin, atom_end,
                                             value, atom_field);
#endif // COLVARS_USE_SOA
  }
  return COLVARS_OK;
}

#endif

#if CMK_SMP && USE_CKLOOP // SMP only

colvarproxy::smp_mode_t colvarproxy_namd::get_smp_mode() const {
  return smp_mode;
}

int colvarproxy_namd::set_smp_mode(smp_mode_t mode) {
  smp_mode = mode;
  return COLVARS_OK;
}


int colvarproxy_namd::smp_loop(int n_items, std::function<int (int)> const &worker)
{
  auto cmkWorker = [&](int start, int end, void * /* result */) {
#if CMK_TRACE_ENABLED
    double before = CmiWallTimer();
#endif
    for (int i = start; i <= end; i++) {
      worker(i);
    }
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(GLOBAL_MASTER_CKLOOP_CALC_ITEM, before, CmiWallTimer());
#endif
  };
  const int numChunks = smp_num_threads() > n_items ?
                        n_items :
                        smp_num_threads();
  cvm::increase_depth();
  CkLoop_Parallelize(numChunks, 0, n_items - 1, cmkWorker, nullptr, CKLOOP_NONE, nullptr);
  cvm::decrease_depth();
  // CkLoop does not support bitwise-OR reduction, so we just return the global error flag
  return cvm::get_error();
}


void calc_cv_biases_smp(int first, int last, void *result, int paramNum, void *param)
{
  colvarproxy_namd *proxy = (colvarproxy_namd *) param;
  colvarmodule *cv = proxy->colvars;
#if CMK_TRACE_ENABLED
  double before = CmiWallTimer();
#endif
  cvm::increase_depth();
  for (int i = first; i <= last; i++) {
    colvarbias *b = (*(cv->biases_active()))[i];
    if (cvm::debug()) {
      cvm::log("["+cvm::to_str(proxy->smp_thread_id())+"/"+cvm::to_str(proxy->smp_num_threads())+
               "]: calc_cv_biases_smp(), first = "+cvm::to_str(first)+
               ", last = "+cvm::to_str(last)+", bias = "+
               b->name+"\n");
    }
    b->update();
  }
  cvm::decrease_depth();
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(GLOBAL_MASTER_CKLOOP_CALC_BIASES,before,CmiWallTimer());
#endif
}


int colvarproxy_namd::smp_biases_loop()
{
  colvarmodule *cv = this->colvars;
  const int numChunks = smp_num_threads() > cv->variables_active_smp()->size() ?
                        cv->variables_active_smp()->size() :
                        smp_num_threads();
  CkLoop_Parallelize(calc_cv_biases_smp, 1, this,
                     numChunks, 0, cv->biases_active()->size()-1);
  return cvm::get_error();
}


void calc_cv_scripted_forces(int paramNum, void *param)
{
  colvarproxy_namd *proxy = (colvarproxy_namd *) param;
  colvarmodule *cv = proxy->colvars;
#if CMK_TRACE_ENABLED
  double before = CmiWallTimer();
#endif
  if (cvm::debug()) {
    cvm::log("["+cvm::to_str(proxy->smp_thread_id())+"/"+cvm::to_str(proxy->smp_num_threads())+
             "]: calc_cv_scripted_forces()\n");
  }
  cv->calc_scripted_forces();
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(GLOBAL_MASTER_CKLOOP_CALC_SCRIPTED_BIASES,before,CmiWallTimer());
#endif
}


int colvarproxy_namd::smp_biases_script_loop()
{
  colvarmodule *cv = this->colvars;
  CkLoop_Parallelize(calc_cv_biases_smp, 1, this,
                     cv->biases_active()->size(), 0, cv->biases_active()->size()-1,
                     1, NULL, CKLOOP_NONE,
                     calc_cv_scripted_forces, 1, this);
  return cvm::get_error();
}

#endif  // #if CMK_SMP && USE_CKLOOP


int colvarproxy_namd::check_replicas_enabled() {
#if CMK_HAS_PARTITION
  return COLVARS_OK;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_namd::replica_index() {
  return CmiMyPartition();
}


int colvarproxy_namd::num_replicas() {
  return CmiNumPartitions();
}


void colvarproxy_namd::replica_comm_barrier() {
  replica_barrier();
}


int colvarproxy_namd::replica_comm_recv(char* msg_data, int buf_len,
                                        int src_rep) {
  DataMessage *recvMsg = NULL;
  replica_recv(&recvMsg, src_rep, CkMyPe());
  CmiAssert(recvMsg != NULL);
  int retval = recvMsg->size;
  if (buf_len >= retval) {
    memcpy(msg_data,recvMsg->data,retval);
  } else {
    retval = 0;
  }
  CmiFree(recvMsg);
  return retval;
}


int colvarproxy_namd::replica_comm_send(char* msg_data, int msg_len,
                                        int dest_rep) {
  replica_send(msg_data, msg_len, dest_rep, CkMyPe());
  return msg_len;
}


/// Request alchemical energy computation every freq steps
int colvarproxy_namd::request_alch_energy_freq(int const freq) {
  // This test is only valid for NAMD3
  if (freq % simparams->computeEnergies) {
    cvm::error("computeEnergies must be a divisor of lambda-dynamics period (" + cvm::to_str(freq) + ").\n");
    return COLVARS_INPUT_ERROR;
  }
  if (!simparams->alchOn) {
    cvm::error("alchOn must be enabled for lambda-dynamics.\n");
    return COLVARS_INPUT_ERROR;
  }
  if (!simparams->alchThermIntOn) {
    cvm::error("alchType must be set to TI for lambda-dynamics.\n");
    return COLVARS_INPUT_ERROR;
  }
  return COLVARS_OK;
}


/// Get value of alchemical lambda parameter from back-end
int colvarproxy_namd::get_alch_lambda(cvm::real* lambda) {
  *lambda = simparams->getCurrentLambda(step);
  return COLVARS_OK;
}


/// Set value of alchemical lambda parameter in back-end
int colvarproxy_namd::send_alch_lambda(void) {
  if (simparams->alchLambdaFreq > 0) {
    cvm::error("Cannot set lambda from Colvars because alchLambdaFreq is enabled. "
                "Either remove biasing forces and extended Lagrangian dynamics on the alchemical coordinate, "
                "or disable alchLambdaFreq.\n");
    return COLVARS_INPUT_ERROR;
  } else {
    simparams->alchLambda = cached_alch_lambda;
  }
  return COLVARS_OK;
}


/// Get energy derivative with respect to lambda
int colvarproxy_namd::get_dE_dlambda(cvm::real* dE_dlambda) {
  // Force data at step zero is garbage in NAMD3, zero in NAMD2
  if (cvm::step_relative() > 0) {
    *dE_dlambda = controller->getTIderivative();
  } else {
    *dE_dlambda = 0.0;
  }
  return COLVARS_OK;
}
