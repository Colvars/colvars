/// -*- c++ -*-

#include <errno.h>

#include "common.h"
#include "fstream_namd.h"
#include "BackEnd.h"
#include "InfoStream.h"
#include "Node.h"
#include "Molecule.h"
#include "PDB.h"
#include "PDBData.h"
#include "ReductionMgr.h"
#include "ScriptTcl.h"

#ifdef NAMD_TCL
#include <tcl.h>
#endif

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "colvarproxy_namd.h"
#include "colvarscript.h"

colvarproxy_namd::colvarproxy_namd()
{
  first_timestep = true;
  system_force_requested = false;
  requestTotalForce(system_force_requested);

  // initialize pointers to NAMD configuration data
  simparams = Node::Object()->simParameters;

  if (cvm::debug())
    iout << "Info: initializing the colvars proxy object.\n" << endi;

  // find the configuration file, if provided
  StringList *config = Node::Object()->configList->find("colvarsConfig");

  // find the input state file
  StringList *input_restart = Node::Object()->configList->find("colvarsInput");
  input_prefix_str = std::string(input_restart ? input_restart->data : "");
  if (input_prefix_str.rfind(".colvars.state") != std::string::npos) {
    // strip the extension, if present
    input_prefix_str.erase(input_prefix_str.rfind(".colvars.state"),
                           std::string(".colvars.state").size());
  }

  // get the thermostat temperature
  if (simparams->rescaleFreq > 0)
    thermostat_temperature = simparams->rescaleTemp;
  else if (simparams->reassignFreq > 0)
    thermostat_temperature = simparams->reassignTemp;
  else if (simparams->langevinOn)
    thermostat_temperature = simparams->langevinTemp;
  else if (simparams->tCoupleOn)
    thermostat_temperature = simparams->tCoupleTemp;
  //else if (simparams->loweAndersenOn)
  //  thermostat_temperature = simparams->loweAndersenTemp;
  else
    thermostat_temperature = 0.0;

  random = Random(simparams->randomSeed);

  // take the output prefixes from the namd input
  output_prefix_str = std::string(simparams->outputFilename);
  restart_output_prefix_str = std::string(simparams->restartFilename);
  restart_frequency_s = simparams->restartFrequency;

  // check if it is possible to save output configuration
  if ((!output_prefix_str.size()) && (!restart_output_prefix_str.size())) {
    fatal_error("Error: neither the final output state file or "
                "the output restart file could be defined, exiting.\n");
  }


#ifdef NAMD_TCL
  have_scripts = true;
  // Store pointer to NAMD's Tcl interpreter
  interp = Node::Object()->getScript()->interp;

  // See is user-scripted forces are defined
  if (Tcl_FindCommand(interp, "calc_colvar_forces", NULL, 0) == NULL) {
    force_script_defined = false;
  } else {
    force_script_defined = true;
  }
#else
  force_script_defined = false;
  have_scripts = false;
#endif


  // initiate module: this object will be the communication proxy
  colvars = new colvarmodule(this);
  cvm::log("Using NAMD interface, version "+
           cvm::to_str(COLVARPROXY_VERSION)+".\n");

  if (config) {
    colvars->read_config_file(config->data);
  }
  colvars->setup_input();
  colvars->setup_output();

  // save to Node for Tcl script access
  Node::Object()->colvars = colvars;


#ifdef NAMD_TCL
  // Construct instance of colvars scripting interface
  script = new colvarscript(this);
#endif


  if (simparams->firstTimestep != 0) {
    cvm::log("Initializing step number as firstTimestep.\n");
    colvars->it = colvars->it_restart = simparams->firstTimestep;
  }

  if (cvm::debug()) {
    cvm::log("colvars_atoms = "+cvm::to_str(colvars_atoms)+"\n");
    cvm::log("colvars_atoms_ncopies = "+cvm::to_str(colvars_atoms_ncopies)+"\n");
    cvm::log("positions = "+cvm::to_str(positions)+"\n");
    cvm::log("total_forces = "+cvm::to_str(total_forces)+"\n");
    cvm::log("applied_forces = "+cvm::to_str(applied_forces)+"\n");
    cvm::log(cvm::line_marker);
  }

  // Initialize reduction object to submit restraint energy as MISC
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  if (cvm::debug())
    iout << "Info: done initializing the colvars proxy object.\n" << endi;
}

/*
  void colvarproxy_namd::construct_cvm(char const  *config_filename)
  // TODO This method might need some refinements for delayed initialization
  // eg. accept config string instead of filename, as below
  //void colvarproxy_namd::construct_cvm (std::string const &config)
  {

  // initiate the colvarmodule, this object will be the communication
  // proxy
  colvars = new colvarmodule(config_filename, this);
  // save to Node for Tcl script access
  Node::Object()->colvars = colvars;

  if (simparams->firstTimestep != 0) {
  cvm::log("Initializing step number as firstTimestep.\n");
  colvars->it = colvars->it_restart = simparams->firstTimestep;
  }

  if (cvm::debug()) {
  cvm::log("colvars_atoms = "+cvm::to_str(colvars_atoms)+"\n");
  cvm::log("colvars_atoms_ncopies = "+cvm::to_str(colvars_atoms_ncopies)+"\n");
  cvm::log("positions = "+cvm::to_str(positions)+"\n");
  cvm::log("total_forces = "+cvm::to_str(total_forces)+"\n");
  cvm::log("applied_forces = "+cvm::to_str(applied_forces)+"\n");
  cvm::log(cvm::line_marker);
  }
  }
*/

colvarproxy_namd::~colvarproxy_namd()
{
  delete reduction;
  if (script != NULL) {
    delete script;
    script = NULL;
  }
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
}

// Reimplemented function from GlobalMaster
void colvarproxy_namd::calculate()
{

  if (first_timestep) {
    first_timestep = false;
  } else {
    // Use the time step number inherited from GlobalMaster
    if ( step - previous_NAMD_step == 1 ) {
      colvars->it++;
    }
    // Other cases could mean:
    // - run 0
    // - beginning of a new run statement
    // then the internal counter should not be incremented
  }
  previous_NAMD_step = step;

  if (cvm::debug()) {
    cvm::log(cvm::line_marker+
             "colvarproxy_namd, step no. "+cvm::to_str(colvars->it)+"\n"+
             "Updating internal data.\n");
  }

  // must delete the forces applied at the previous step: they have
  // already been used and copied to other memory locations
  modifyForcedAtoms().resize(0);
  modifyAppliedForces().resize(0);

  // prepare the local arrays to contain the sorted copies of the NAMD
  // arrays
  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    positions[i] = cvm::rvector(0.0, 0.0, 0.0);
    total_forces[i] = cvm::rvector(0.0, 0.0, 0.0);
    applied_forces[i] = cvm::rvector(0.0, 0.0, 0.0);
  }

  // sort the positions array
  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    bool found_position = false;
    AtomIDList::const_iterator a_i = this->getAtomIdBegin();
    AtomIDList::const_iterator a_e = this->getAtomIdEnd();
    PositionList::const_iterator p_i = this->getAtomPositionBegin();
    for ( ; a_i != a_e; ++a_i, ++p_i ) {
      if ( *a_i == colvars_atoms[i] ) {
        found_position = true;
        Position const &namd_pos = *p_i;
        positions[i] = cvm::rvector(namd_pos.x, namd_pos.y, namd_pos.z);
        break;
      }
    }
    if (!found_position)
      cvm::fatal_error("Error: cannot find the position of atom "+
                       cvm::to_str(colvars_atoms[i]+1)+"\n");
  }


  if (system_force_requested && cvm::step_relative() > 0) {

    // sort the array of total forces from the previous step (but only
    // do it if there *is* a previous step!)
    for (size_t i = 0; i < colvars_atoms.size(); i++) {
      bool found_total_force = false;
      //found_total_force = false;
      AtomIDList::const_iterator a_i = this->getForceIdBegin();
      AtomIDList::const_iterator a_e = this->getForceIdEnd();
      PositionList::const_iterator f_i = this->getTotalForce();
      for ( ; a_i != a_e; ++a_i, ++f_i ) {
        if ( *a_i == colvars_atoms[i] ) {
          found_total_force = true;
          Vector const &namd_force = *f_i;
          total_forces[i] = cvm::rvector(namd_force.x, namd_force.y, namd_force.z);
          //           if (cvm::debug())
          //             cvm::log ("Found the total force of atom "+
          //                       cvm::to_str (colvars_atoms[i]+1)+", which is "+
          //                       cvm::to_str (total_forces[i])+".\n");
          break;
        }
      }
      if (!found_total_force)
        cvm::fatal_error("Error: system forces were requested, but total force on atom "+
                         cvm::to_str(colvars_atoms[i]+1) + " was not\n"
                         "found. The most probable cause is combination of energy minimization with a\n"
                         "biasing method that requires MD (e.g. ABF). Always run minimization\n"
                         "and ABF separately.");
    }

    // do the same for applied forces
    for (size_t i = 0; i < colvars_atoms.size(); i++) {
      AtomIDList::const_iterator a_i = this->getLastAtomsForcedBegin();
      AtomIDList::const_iterator a_e = this->getLastAtomsForcedEnd();
      PositionList::const_iterator f_i = this->getLastForcesBegin();
      for ( ; a_i != a_e; ++a_i, ++f_i ) {
        if ( *a_i == colvars_atoms[i] ) {
          Vector const &namd_force = *f_i;
          if (cvm::debug())
            cvm::log("Found a force applied to atom "+
                     cvm::to_str(colvars_atoms[i]+1)+": "+
                     cvm::to_str(cvm::rvector(namd_force.x, namd_force.y, namd_force.z))+
                     "; current total is "+
                     cvm::to_str(applied_forces[i])+".\n");
          applied_forces[i] += cvm::rvector(namd_force.x, namd_force.y, namd_force.z);
        }
      }
    }
  }

  // call the collective variable module
  if (colvars->calc() != COLVARS_OK) {
    fatal_error("");
  }
  // send MISC energy
  reduction->submit();

  // NAMD does not destruct GlobalMaster objects, so we must remember
  // to write all output files at the end of the run
  if (step == simparams->N) {
    colvars->write_output_files();
  }
}

// Callback functions

int colvarproxy_namd::run_force_callback() {
#ifdef NAMD_TCL
  std::string cmd = std::string("calc_colvar_forces ")
    + cvm::to_str(cvm::step_absolute());
  int err = Tcl_Eval(interp, cmd.c_str());
  if (err != TCL_OK) {
    cvm::log(std::string("Error while executing calc_colvar_forces:\n"));
    cvm::error(Tcl_GetStringResult(interp));
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}

int colvarproxy_namd::run_colvar_callback(std::string const &name,
                                          std::vector<const colvarvalue *> const &cvc_values,
                                          colvarvalue &value)
{
#ifdef NAMD_TCL
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
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}

int colvarproxy_namd::run_colvar_gradient_callback(std::string const &name,
                                                   std::vector<const colvarvalue *> const &cvc_values,
                                                   std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
#ifdef NAMD_TCL
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
    cvm::error("Error parsing list of gradient values from script: found "
        + cvm::to_str(n) + " values instead of " + cvm::to_str(gradient.size()));
    return COLVARS_ERROR;
  }
  for (i = 0; i < gradient.size(); i++) {
    std::istringstream is(Tcl_GetString(list[i]));
    if (gradient[i].from_simple_string(is.str()) != COLVARS_OK) {
      cvm::log("Gradient matrix size: " + cvm::to_str(gradient[i].size()));
      cvm::log("Gradient string: " + cvm::to_str(Tcl_GetString(list[i])));
      cvm::error("Error parsing gradient value from script");
      return COLVARS_ERROR;
    }
  }
  return (err == TCL_OK) ? COLVARS_OK : COLVARS_ERROR;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


void colvarproxy_namd::add_energy(cvm::real energy)
{
  reduction->item(REDUCTION_MISC_ENERGY) += energy;
}

void colvarproxy_namd::request_system_force(bool yesno)
{
  system_force_requested = yesno;
  requestTotalForce(system_force_requested);
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
  // In NAMD, all errors are fatal
  fatal_error(message);
}


void colvarproxy_namd::fatal_error(std::string const &message)
{
  log(message);
  if (errno) log(strerror(errno));
  if (!cvm::debug())
    log("If this error message is unclear, "
        "try recompiling with -DCOLVARS_DEBUG.\n");
  if (errno) {
    NAMD_err("Error in the collective variables module");
  } else {
    NAMD_die("Error in the collective variables module: exiting.\n");
  }
}


void colvarproxy_namd::exit(std::string const &message)
{
  cvm::log(message);
  BackEnd::exit();
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
               pdb_field_str+"\".\n", INPUT_ERROR);
  }

  return pdb_field;
}


int colvarproxy_namd::load_coords(char const *pdb_filename,
                                  std::vector<cvm::atom_pos> &pos,
                                  const std::vector<int> &indices,
                                  std::string const &pdb_field_str,
                                  double const pdb_field_value)
{
  if (pdb_field_str.size() == 0 && indices.size() == 0) {
    cvm::fatal_error("Bug alert: either PDB field should be defined or list of "
                     "atom IDs should be available when loading atom coordinates!\n");
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
        if (ipdb != *current_index) {
          // Skip atoms not in the list
          continue;
        } else {
          current_index++;
        }
      }

      if (!pos_allocated) {
        pos.push_back(cvm::atom_pos(0.0, 0.0, 0.0));
      } else if (ipos >= pos.size()) {
        cvm::fatal_error("Error: the PDB file \""+
                         std::string(pdb_filename)+
                         "\" contains coordinates for "
                         "more atoms than needed.\n");
      }

      pos[ipos] = cvm::atom_pos((pdb->atom(ipdb))->xcoor(),
                                (pdb->atom(ipdb))->ycoor(),
                                (pdb->atom(ipdb))->zcoor());
      ipos++;
      if (!use_pdb_field && current_index == indices.end())
        break;
    }

    if ((ipos < pos.size()) || (current_index != indices.end()))
      cvm::fatal_error("Error: the number of records in the PDB file \""+
                       std::string(pdb_filename)+
                       "\" does not appear to match either the total number of atoms,"+
                       " or the number of coordinates requested at this point("+
                       cvm::to_str(pos.size())+").\n");

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


int colvarproxy_namd::load_atoms(char const *pdb_filename,
                                 std::vector<cvm::atom> &atoms,
                                 std::string const &pdb_field_str,
                                 double const pdb_field_value)
{
  if (pdb_field_str.size() == 0)
    cvm::fatal_error("Error: must define which PDB field to use "
                     "in order to define atoms from a PDB file.\n");

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

    atoms.push_back(cvm::atom(ipdb+1));
  }

  delete pdb;
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


std::ostream * colvarproxy_namd::output_stream(std::string const &output_name)
{
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osni == output_name) {
      return *osi;
    }
  }
  output_stream_names.push_back(output_name);
  this->backup_file(output_name.c_str());
  ofstream_namd * os = new ofstream_namd(output_name.c_str());
  if (!os->is_open()) {
    cvm::error("Error: cannot write to file \""+output_name+"\".\n",
               FILE_ERROR);
  }
  output_files.push_back(os);
  return os;
}

int colvarproxy_namd::close_output_stream(std::string const &output_name)
{
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osni == output_name) {
      ((ofstream_namd *) *osi)->close();
      output_files.erase(osi);
      output_stream_names.erase(osni);
      return COLVARS_OK;
    }
  }
  return COLVARS_ERROR;
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


size_t colvarproxy_namd::init_namd_atom(AtomID const &aid)
{
  modifyRequestedAtoms().add(aid);
  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    if (colvars_atoms[i] == aid) {
      // this atom id was already recorded
      colvars_atoms_ncopies[i] += 1;
      return i;
    }
  }

  // allocate a new slot for this atom
  colvars_atoms_ncopies.push_back(1);
  colvars_atoms.push_back(aid);
  positions.push_back(cvm::rvector());
  total_forces.push_back(cvm::rvector());
  applied_forces.push_back(cvm::rvector());

  return (colvars_atoms.size()-1);
}

// atom member functions, NAMD specific implementations

cvm::atom::atom(int const &atom_number)
{
  // NAMD internal numbering starts from zero
  AtomID const aid(atom_number-1);

  if (cvm::debug())
    cvm::log("Adding atom "+cvm::to_str(aid+1)+
             " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= Node::Object()->molecule->numAtoms) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n");
    return;
  }
  this->index = ((colvarproxy_namd *) cvm::proxy)->init_namd_atom(aid);
  if (cvm::debug())
    cvm::log("The index of this atom in the colvarproxy_namd arrays is "+
             cvm::to_str(this->index)+".\n");
  this->id = aid;
  this->mass = Node::Object()->molecule->atommass(aid);
  this->reset_data();
}


/// For AMBER topologies, the segment id is automatically set to
/// "MAIN" (the segment id assigned by NAMD's AMBER topology parser),
/// and is therefore optional when an AMBER topology is used
cvm::atom::atom(cvm::residue_id const &residue,
                std::string const     &atom_name,
                std::string const     &segment_id)
{
  AtomID const aid =
    (segment_id.size() ?
     Node::Object()->molecule->get_atom_from_name(segment_id.c_str(),
                                                  residue,
                                                  atom_name.c_str()) :
     Node::Object()->molecule->get_atom_from_name("MAIN",
                                                  residue,
                                                  atom_name.c_str()));


  if (cvm::debug())
    cvm::log("Adding atom \""+
             atom_name+"\" in residue "+
             cvm::to_str(residue)+
             " (index "+cvm::to_str(aid)+
             ") for collective variables calculation.\n");

  if (aid < 0) {
    // get_atom_from_name() has returned an error value
    cvm::fatal_error("Error: could not find atom \""+
                     atom_name+"\" in residue "+
                     cvm::to_str(residue)+
                     ( (segment_id != "MAIN") ?
                       (", segment \""+segment_id+"\"") :
                       ("") )+
                     "\n");
  }

  this->index = ((colvarproxy_namd *) cvm::proxy)->init_namd_atom(aid);
  if (cvm::debug())
    cvm::log("The index of this atom in the colvarproxy_namd arrays is "+
             cvm::to_str(this->index)+".\n");
  this->id = aid;
  this->mass = Node::Object()->molecule->atommass(aid);
  this->reset_data();
}


// copy constructor
cvm::atom::atom(cvm::atom const &a)
  : index(a.index), id(a.id), mass(a.mass)
{
  // init_namd_atom() has already been called by a's constructor, no
  // need to call it again

  // need to increment the counter anyway
  colvarproxy_namd *gm = (colvarproxy_namd *) cvm::proxy;
  gm->colvars_atoms_ncopies[this->index] += 1;
}


cvm::atom::~atom()
{
  if (this->index >= 0) {
    colvarproxy_namd *gm = (colvarproxy_namd *) cvm::proxy;
    if (gm->colvars_atoms_ncopies[this->index] > 0)
      gm->colvars_atoms_ncopies[this->index] -= 1;
  }
}


void cvm::atom::read_position()
{
  colvarproxy_namd const * const gm = (colvarproxy_namd *) cvm::proxy;
  this->pos = gm->positions[this->index];
}


void cvm::atom::read_velocity()
{
  cvm::fatal_error("Error: NAMD does not have yet a way to communicate "
                   "atom velocities to the colvars.\n");
}


void cvm::atom::read_system_force()
{
  colvarproxy_namd const * const gm = (colvarproxy_namd *) cvm::proxy;
  this->system_force = gm->total_forces[this->index] - gm->applied_forces[this->index];
}


void cvm::atom::apply_force(cvm::rvector const &new_force)
{
  colvarproxy_namd *gm = (colvarproxy_namd *) cvm::proxy;
  gm->modifyForcedAtoms().add(this->id);
  gm->modifyAppliedForces().add(Vector(new_force.x, new_force.y, new_force.z));
}

