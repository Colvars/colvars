/// -*- c++ -*-
/* Jeff Comer's tests to see if he can link GROMACS and Colvars */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cerrno>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//#include "gromacs/fileio/futil.h"
//#include "index.h"
//#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vec.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/network.h"
//#include "gromacs/fileio/filenm.h"
#include <string.h>
#include "gromacs/utility/smalloc.h"
//#include "pull.h"
//#include "xvgr.h"
//#include "names.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/gmx_ga2la.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "colvars_potential.h"
#include "colvarproxy_gromacs.h"
#include "gromacs/random/random.h"
#include "gromacs/utility/futil.h"

// Jeff Comer's tests to see if he can link GROMACS and Colvars
#include "colvars_potential.h"
#include "colvarproxy_gromacs.h"

// The global colvars proxy object, awkwardly initialized on
// the first call to colvars_potential.
colvarproxy_gromacs colvars_global_proxy;
bool colvars_global_is_first = true;


real colvars_potential(t_inputrec *gmx_inp, t_mdatoms *md, t_pbc *pbc,
		       gmx_int64_t step, rvec *x, rvec *f, tensor vir) {

  // Update some things.
  // Get the current periodic boundary conditions.
  colvars_global_proxy.gmx_pbc = (*pbc);
  colvars_global_proxy.gmx_atoms = md;

  // Get the thermostat temperature.
  // NOTE: Considers only the first temperature coupling group!
  // I'm not sure if this can change during the simulation, so
  // putting it every step to be safe.
  if (gmx_inp->opts.ref_t[0] > 0)
    colvars_global_proxy.set_temper(gmx_inp->opts.ref_t[0]);
  else
    colvars_global_proxy.set_temper(gmx_inp->opts.ref_t[0]); // FIXME duplicate or 'if' case above

  // Initialize if this is the first call.
  if (colvars_global_is_first) {
    colvars_global_proxy.init(gmx_inp, step);
    colvars_global_is_first = false;
  }

  // colvars computation
  return colvars_global_proxy.calculate(step, x, f, vir);
}

// Taken from colvarproxy_lammps.
// safely move filename to filename.extension
static int my_backup_file(const char *filename, const char *extension)
{
  struct stat sbuf;
  if (stat(filename, &sbuf) == 0) {
    if (!extension) extension = ".BAK";
    char *backup = new char[strlen(filename)+strlen(extension)+1];
    strcpy(backup, filename);
    strcat(backup, extension);
#if defined(_WIN32) && !defined(__CYGWIN__)
    remove(backup);
#endif
    if (rename(filename,backup)) {
      char *sys_err_msg = strerror(errno);
      if (!sys_err_msg)  sys_err_msg = (char *) "(unknown error)";
      fprintf(stderr,"Error renaming file %s to %s: %s\n",
              filename, backup, sys_err_msg);
      delete [] backup;
      return COLVARS_ERROR;
    }
    delete [] backup;
  }
  return COLVARS_OK;
}



//************************************************************
// colvarproxy_gromacs
colvarproxy_gromacs::colvarproxy_gromacs() : colvarproxy() {
}

// Colvars Initialization
void colvarproxy_gromacs::init(t_inputrec *gmx_inp, gmx_int64_t step) {
  // Initialize colvars.
  first_timestep = true;
  total_force_requested = false;

  // User-scripted forces are not available in GROMACS
  force_script_defined = false;
  have_scripts = false;

  angstrom_value = 0.1;

  // GROMACS random number generation.
  // Seed with the mdp parameter ld_seed, the Langevin dynamics seed.
  rando = gmx_rng_init(gmx_inp->ld_seed);

  // For expediency, we are using a kludgy input/output structure
  //
  // gmx_inp->userint1 is the colvars state:
  // 0, no colvars
  // 1, colvars without restart
  // 2, colvars with restart (from colvars.state.restart)
  //
  // gmx_inp->userint2 is the config file index
  // the colvars config file is assumed to be called colvars${userint2}.colvars
  // the output files will have the prefix colvars${userint2}
  //
  // gmx_inp->userint3 is the input file index
  // the input state file is assumed to have the prefix colvars${userint3}

  switch(gmx_inp->userint1) {
  case 0:
    gmx_fatal(FARGS,"colvars_potential called for userint1=0\n");
    break;
  case 1:
    colvars_restart = false;
    break;
  case 2:
    colvars_restart = true;
    break;
  default:
    gmx_fatal(FARGS,"userint1 must be 1 (colvars on) or 2 (restarted colvars), not %d\n", gmx_inp->userint1);
    break;
  }
  // Right now the input/output prefixes almost hardcoded (expect for an integer)
  char prefix[256];
  snprintf(prefix,256,"colvars%d",gmx_inp->userint2);
  config_file = std::string(prefix).append(".colvars");
  // Output
  output_prefix_str = std::string(prefix);
  restart_output_prefix_str = std::string(prefix).append(".restart");
  // State file
  if (colvars_restart) {
    char input_state_prefix[256];
    snprintf(input_state_prefix,256,"colvars%d",gmx_inp->userint3);
    input_prefix_str = std::string(input_state_prefix);
  }

  // Get some parameters from GROMACS
  restart_frequency_s = gmx_inp->nstxout;
  timestep = gmx_inp->delta_t;

  // initiate module: this object will be the communication proxy
  colvars = new colvarmodule (this);
  cvm::log("Using GROMACS interface, version "+
	   cvm::to_str(COLVARPROXY_VERSION)+".\n");
  colvars->read_config_file(config_file.c_str());
  colvars->setup_input();
  colvars->setup_output();

  if (step != 0) {
    cvm::log("Initializing step number to "+cvm::to_str(step)+".\n");
    colvars->it = colvars->it_restart = step;
  }

  if (cvm::debug()) {
    cvm::log ("colvars_atoms = "+cvm::to_str (colvars_atoms)+"\n");
    cvm::log ("colvars_atoms_ncopies = "+cvm::to_str (colvars_atoms_ncopies)+"\n");
    cvm::log ("positions = "+cvm::to_str (positions)+"\n");
    cvm::log ("total_forces = "+cvm::to_str (total_forces)+"\n");
    cvm::log ("applied_forces = "+cvm::to_str (applied_forces)+"\n");
    cvm::log (cvm::line_marker);
  }

  if (cvm::debug())
    log("done initializing the colvars proxy object.\n");
} // End colvars initialization.

// We don't really know when GROMACS will die, but
// since the object has file scope, this should get called.
colvarproxy_gromacs::~colvarproxy_gromacs()
{
  if (colvars != NULL) {
    colvars->write_output_files();
    delete colvars;
    colvars = NULL;
  }
}

void colvarproxy_gromacs::set_temper(double temper) {
  thermostat_temperature = temper;
}

// GROMACS uses nanometers and kJ/mol internally
cvm::real colvarproxy_gromacs::backend_angstrom_value() { return 0.1; }

// From Gnu units
// $ units -ts 'k' 'kJ/mol/K/avogadro'
// 0.0083144599 with v2.16 (older value 0.0083144621)
cvm::real colvarproxy_gromacs::boltzmann() { return 0.0083144599; }

// Temperature of the simulation (K)
cvm::real colvarproxy_gromacs::temperature() {
  return thermostat_temperature;
}

// Time step of the simulation (fs)
// GROMACS uses picoseconds.
cvm::real colvarproxy_gromacs::dt() { return 1000.0*timestep; }

cvm::real colvarproxy_gromacs::rand_gaussian() {
  return gmx_rng_gaussian_real(rando);
}

void colvarproxy_gromacs::request_total_force (bool yesno) {
  total_force_requested = yesno;
}

size_t colvarproxy_gromacs::restart_frequency() {
  return restart_frequency_s;
}

// **************** PERIODIC BOUNDARY CONDITIONS ****************
//  Get the PBC-aware distance vector between two positions
cvm::rvector colvarproxy_gromacs::position_distance (cvm::atom_pos const &pos1,
				cvm::atom_pos const &pos2) {
  rvec r1, r2, dr;
  r1[0] = pos1.x;
  r1[1] = pos1.y;
  r1[2] = pos1.z;
  r2[0] = pos2.x;
  r2[1] = pos2.y;
  r2[2] = pos2.z;

  pbc_dx(&gmx_pbc, r1, r2, dr);
  return cvm::atom_pos( dr[0], dr[1], dr[2] );
}

// The position to look for the closest periodic image
// ref_pos The reference position
void colvarproxy_gromacs::select_closest_image (cvm::atom_pos &pos,
			   cvm::atom_pos const &ref_pos) {
  rvec r1, r2, dr;
  r1[0] = pos.x;
  r1[1] = pos.y;
  r1[2] = pos.z;
  r2[0] = ref_pos.x;
  r2[1] = ref_pos.y;
  r2[2] = ref_pos.z;
  pbc_dx_aiuc(&gmx_pbc, r1, r2, dr);

  // dr is the closest distance vector.
  pos.x = r2[0]+dr[0];
  pos.y = r1[1]+dr[1];
  pos.z = r1[2]+dr[2];
}

void colvarproxy_gromacs::log (std::string const &message) {
  printf("colvars: %s", message.c_str());
}

void colvarproxy_gromacs::error (std::string const &message) {
  // In GROMACS, all errors are fatal.
  fatal_error (message);
}

void colvarproxy_gromacs::fatal_error (std::string const &message) {
  log(message);
  if (!cvm::debug())
    log("If this error message is unclear, "
	"try recompiling with -DCOLVARS_DEBUG.\n");
  gmx_fatal(FARGS,"Error in collective variables module.\n");
}

void colvarproxy_gromacs::exit (std::string const &message) {
  gmx_fatal(FARGS,"SUCCESS: %s\n", message.c_str());
}

int colvarproxy_gromacs::load_atoms (char const *filename,
		std::vector<cvm::atom> &atoms,
		std::string const &pdb_field,
		double const pdb_field_value) {
  cvm::error("Selecting collective variable atoms "
		   "from a PDB file is currently not supported.\n");
  return COLVARS_NOT_IMPLEMENTED;
}

int colvarproxy_gromacs::load_coords (char const *filename,
                                    std::vector<cvm::atom_pos> &pos,
                                    const std::vector<int> &indices,
                                    std::string const &pdb_field_str,
                                    double const pdb_field_value) {
  cvm::error("Selecting collective variable atoms "
		   "from a PDB file is currently not supported.\n");
  return COLVARS_NOT_IMPLEMENTED;
}

int colvarproxy_gromacs::set_unit_system(std::string const &units_in, bool /*check_only*/)
{
  if (units_in != "gromacs") {
    cvm::error("Specified unit system \"" + units_in + "\" is unsupported in Gromacs. Supported units are \"gromacs\" (nm, kJ/mol).\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}

int colvarproxy_gromacs::backup_file (char const *filename)
{
  if (std::string(filename).rfind(std::string(".colvars.state"))
      != std::string::npos) {
    return my_backup_file(filename, ".old");
  } else {
    // GROMACS has its own way to avoid overwriting files.
    //if (make_backup(filename)) return COLVARS_OK;
    //else return FILE_ERROR;
    make_backup(filename);
  }
  return COLVARS_OK;
}


// trigger colvars computation
double colvarproxy_gromacs::calculate(gmx_int64_t step, const rvec *x, rvec *f, tensor vir) {
  if (first_timestep) {
    first_timestep = false;
  } else {
    // Use the time step number inherited from GROMACS
    if ( step - previous_gmx_step == 1 )
      colvars->it++;
    // Other cases?
  }
  previous_gmx_step = step;

  if (cvm::debug()) {
    cvm::log(cvm::line_marker+
             "colvarproxy_gromacs, step no. "+cvm::to_str(colvars->it)+"\n"+
             "Updating internal data.\n");
  }

  // backup applied forces if necessary to calculate total forces
  //if (total_force_requested)
  //  previous_applied_forces = applied_forces;

  // Zero the forces on the atoms, so that they can be accumulated by the colvars.
  for (size_t i = 0; i < applied_forces.size(); i++) {
    applied_forces[i].x = applied_forces[i].y = applied_forces[i].z = 0.0;
  }

  // Get the atom positions from the Gromacs array.
  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    size_t aid = colvars_atoms[i];
    if (aid >= gmx_atoms->nr) {
      cvm::fatal_error("Error: Atom index "+cvm::to_str(aid)+" not found in GROMACS data structure containing "+cvm::to_str(gmx_atoms->nr)+" atoms");
    }
    positions[i] = cvm::rvector(x[aid][0], x[aid][1], x[aid][2]);
  }

  // Get total forces if required.
  if (total_force_requested && cvm::step_relative() > 0) {
     for (size_t i = 0; i < colvars_atoms.size(); i++) {
       size_t aid = colvars_atoms[i];
       // We already checked above that gmx_atoms->nr < aid.
       // The change of sign appears necessary.
       total_forces[i] = cvm::rvector(-f[aid][0], -f[aid][1], -f[aid][2]);
     }
  }

  bias_energy = 0.0;
  // Call the collective variable module to fill applied_forces
  if (colvars->calc() != COLVARS_OK) {
    cvm::fatal_error("");
  }

  // Pass the applied forces back to GROMACS.
  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    size_t aid = colvars_atoms[i];
    // We already checked above that gmx_atoms->nr < aid.
    f[aid][0] -= applied_forces[i].x;
    f[aid][1] -= applied_forces[i].y;
    f[aid][2] -= applied_forces[i].z;
  }

  // We should probably update the virial.

  return bias_energy;
}


// Pass restraint energy value for current timestep to MD engine
void colvarproxy_gromacs::add_energy (cvm::real energy)
{
  bias_energy += energy;
}

// **************** ATOMS ****************
// Most of the following was patterned after colvarproxy_lammps.

int colvarproxy_gromacs::init_gromacs_atom(const int &aid, cvm::atom *atom)
{
  atom->id = aid;
  atom->mass = gmx_atoms->massT[aid];

  if ( (aid < 0) || (aid >= gmx_atoms->nr) ) {
    cvm::fatal_error ("Error: invalid atom number specified, "+
                      cvm::to_str (aid+1)+"\n");
  }

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

  return colvars_atoms.size()-1;
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

cvm::atom::atom(const int &atom_number)
{
  // GROMACS uses zero-based arrays.
  int aid = atom_number-1;

  if (cvm::debug())
    cvm::log ("Adding atom "+cvm::to_str (aid+1)+
              " for collective variables calculation.\n");

  this->index = ((colvarproxy_gromacs *) cvm::proxy)->init_gromacs_atom (aid,this);
  if (cvm::debug())
    cvm::log ("The index of this atom in the colvarproxy_gromacs arrays is "+
              cvm::to_str (this->index)+".\n");
  this->reset_data();
}

/// For AMBER topologies, the segment id is automatically set to
/// "MAIN" (the segment id assigned by NAMD's AMBER topology parser),
/// and is therefore optional when an AMBER topology is used
cvm::atom::atom(cvm::residue_id const &residue,
                std::string const     &atom_name,
                std::string const     &segment_id)
{
  cvm::fatal_error("Error: Creating collective variable atoms "
                   "from a PDB file is currently not supported.\n");
}

// copy constructor
cvm::atom::atom(cvm::atom const &a)
  : index(a.index), id(a.id), mass(a.mass)
{
  // init_gromacs_atom() has already been called by a's constructor, no
  // need to call it again

  // need to increment the counter anyway
  colvarproxy_gromacs *cp = (colvarproxy_gromacs *) cvm::proxy;
  cp->colvars_atoms_ncopies[this->index] += 1;
}

cvm::atom::~atom()
{
  if (this->index >= 0) {
    colvarproxy_gromacs *cp = (colvarproxy_gromacs *) cvm::proxy;
    if (cp->colvars_atoms_ncopies[this->index] > 0)
      cp->colvars_atoms_ncopies[this->index] -= 1;
  }
}

void cvm::atom::read_position()
{
  colvarproxy_gromacs const * const cp = (colvarproxy_gromacs *) cvm::proxy;
  this->pos = cp->positions[this->index];
}

void cvm::atom::read_velocity()
{
  cvm::fatal_error("Error: read_velocity is not yet implemented.\n");
}

void cvm::atom::read_total_force()
{
  colvarproxy_gromacs const * const cp = (colvarproxy_gromacs *) cvm::proxy;
  this->total_force.x = cp->total_forces[this->index].x;
  this->total_force.y = cp->total_forces[this->index].y;
  this->total_force.z = cp->total_forces[this->index].z;
}

void cvm::atom::apply_force(cvm::rvector const &new_force)
{
  colvarproxy_gromacs *cp = (colvarproxy_gromacs *) cvm::proxy;
  cp->applied_forces[this->index].x += new_force.x;
  cp->applied_forces[this->index].y += new_force.y;
  cp->applied_forces[this->index].z += new_force.z;
}

