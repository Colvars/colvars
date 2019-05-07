/// -*- c++ -*-

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cerrno>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"

#include "colvarproxy_gromacs.h"


//************************************************************
// colvarproxy_gromacs
colvarproxy_gromacs::colvarproxy_gromacs() : colvarproxy() {}

// Colvars Initialization
void colvarproxy_gromacs::init(t_inputrec *ir, int64_t step,t_mdatoms *md,
                               const std::string &prefix,
                               gmx::ArrayRef<const std::string> filenames_config,
                               const std::string &filename_restart) {

  version_int = get_version_from_string(COLVARPROXY_VERSION);

  if (cvm::debug())
    log("Initializing the colvars proxy object.\n");

  // Initialize colvars.
  first_timestep = true;
  total_force_requested = false;
  restart_frequency_s = 0;

  // User-scripted forces are not available in GROMACS
  force_script_defined = false;
  have_scripts = false;

  // Get the thermostat temperature.
  // NOTE: Considers only the first temperature coupling group!
  thermostat_temperature = ir->opts.ref_t[0];

  // GROMACS random number generation.
  // Seed with the mdp parameter ld_seed, the Langevin dynamics seed.
  rng.seed(ir->ld_seed);

  /// Handle input filenames and prefix/suffix for colvars files.
  ///
  /// filename_config is the colvars configuration file collected from "-colvars" option.
  /// The output prefix will be the prefix of Gromacs log filename.
  /// or "output" otherwise.
  ///
  /// For restart, 'filename_restart' is the colvars input file for restart,
  /// set by the "-cv_restart" option. It will be NULL otherwise.
  ///

  if(!prefix.empty())
  {
    output_prefix_str = prefix;
  }
  else {
    output_prefix_str = "output";
  }

  restart_output_prefix_str = prefix + ".restart";

  colvars_restart = false;

  if(!filename_restart.empty())
  {
    colvars_restart = true;
    input_prefix_str = filename_restart;
  }

  // Retrieve masses and charges from input file
  updated_masses_ = updated_charges_ = true;

  // Get some parameters from GROMACS
  timestep = ir->delta_t;
  gmx_atoms = md;

  // initiate module: this object will be the communication proxy
  colvars = new colvarmodule (this);
  cvm::log("Using GROMACS interface, version "+
	   cvm::to_str(COLVARPROXY_VERSION)+".\n");

  auto i = filenames_config.begin();
  for(; i != filenames_config.end(); ++i) {
      colvars->read_config_file(i->c_str());
  }

  colvars->setup();
  colvars->setup_input();
  colvars->setup_output();

  if (step != 0) {
    cvm::log("Initializing step number to "+cvm::to_str(step)+".\n");
    colvars->it = colvars->it_restart = step;
  }

  if (cvm::debug()) {
    cvm::log ("atoms_ids = "+cvm::to_str (atoms_ids)+"\n");
    cvm::log ("atoms_ncopies = "+cvm::to_str (atoms_ncopies)+"\n");
    cvm::log ("positions = "+cvm::to_str (atoms_positions)+"\n");
    cvm::log ("total_forces = "+cvm::to_str (atoms_total_forces)+"\n");
    cvm::log ("atoms_new_colvar_forces = "+cvm::to_str (atoms_new_colvar_forces)+"\n");
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
    colvars->write_restart_file(output_prefix_str+".colvars.state");
    colvars->write_output_files();
    delete colvars;
    colvars = NULL;
  }
}

void colvarproxy_gromacs::set_temper(double temper)
{
  thermostat_temperature = temper;
}

// GROMACS uses nanometers.
cvm::real colvarproxy_gromacs::unit_angstrom() { return 0.1; }

// From Gnu units
// $ units -ts 'k' 'kJ/mol/K/avogadro'
// 0.0083144621
cvm::real colvarproxy_gromacs::boltzmann() { return 0.0083144621; }

// Temperature of the simulation (K)
cvm::real colvarproxy_gromacs::temperature()
{
  return thermostat_temperature;
}

// Time step of the simulation (fs)
// GROMACS uses picoseconds.
cvm::real colvarproxy_gromacs::dt() { return 1000.0*timestep; }

cvm::real colvarproxy_gromacs::rand_gaussian()
{
  return  normal_distribution(rng);
}

void colvarproxy_gromacs::request_total_force (bool yesno)
{
  total_force_requested = yesno;
}

size_t colvarproxy_gromacs::restart_frequency()
{
  return restart_frequency_s;
}

// **************** PERIODIC BOUNDARY CONDITIONS ****************
//  Get the PBC-aware distance vector between two positions
cvm::rvector colvarproxy_gromacs::position_distance (cvm::atom_pos const &pos1,
                                                     cvm::atom_pos const &pos2) const
{
  rvec r1, r2, dr;
  r1[0] = pos1.x;
  r1[1] = pos1.y;
  r1[2] = pos1.z;
  r2[0] = pos2.x;
  r2[1] = pos2.y;
  r2[2] = pos2.z;

  pbc_dx(&gmx_pbc, r2, r1, dr);
  return cvm::atom_pos( dr[0], dr[1], dr[2] );
}


void colvarproxy_gromacs::log (std::string const &message)
{
  // Gromacs prints messages on the stderr FILE.
  fprintf(stderr, "colvars: %s", message.c_str());
}

void colvarproxy_gromacs::error (std::string const &message)
{
  // In GROMACS, all errors are fatal.
  fatal_error (message);
}

void colvarproxy_gromacs::fatal_error (std::string const &message)
{
  log(message);
  if (!cvm::debug())
    log("If this error message is unclear, "
	"try recompiling with -DCOLVARS_DEBUG.\n");
  gmx_fatal(FARGS,"Error in collective variables module.\n");
}

void colvarproxy_gromacs::exit (std::string const &message)
{
  gmx_fatal(FARGS,"SUCCESS: %s\n", message.c_str());
}

int colvarproxy_gromacs::load_atoms (char const *filename, std::vector<cvm::atom> &atoms,
                                     std::string const &pdb_field, double const pdb_field_value)
{
  cvm::error("Selecting collective variable atoms "
		   "from a PDB file is currently not supported.\n");
  return COLVARS_NOT_IMPLEMENTED;
}

int colvarproxy_gromacs::load_coords (char const *filename, std::vector<cvm::atom_pos> &pos,
                                      const std::vector<int> &indices, std::string const &pdb_field_str,
                                      double const pdb_field_value)
{
  cvm::error("Selecting collective variable atoms "
		   "from a PDB file is currently not supported.\n");
  return COLVARS_NOT_IMPLEMENTED;
}

int colvarproxy_gromacs::backup_file (char const *filename)
{
  // Incremental gromacs backup system will be use only for those file
  if (std::string(filename).rfind(std::string(".colvars.traj")) != std::string::npos) {

    // GROMACS function
    make_backup(filename);

  // Otherwise, just keep one backup.
  } else {

    //Handle filename of the backup file
    const char *extension = ".old";
    char *backup = new char[strlen(filename)+strlen(extension)+1];
    strcpy(backup, filename);
    strcat(backup, extension);

    gmx_file_copy(filename, backup, FALSE);

    delete [] backup;

  }
  return COLVARS_OK;
}

real colvarproxy_gromacs::colvars_potential(const t_mdatoms *md, t_pbc *pbc,
		                   int64_t step, gmx::ArrayRef<const gmx::RVec> x, gmx::ForceWithVirial *force)
{
  // Update some things.
  // Get the current periodic boundary conditions.
  gmx_pbc = (*pbc);
  gmx_atoms = md;

  // colvars computation
  return calculate(step, x, force);
}

// trigger colvars computation
// TODO: compute the virial contribution
double colvarproxy_gromacs::calculate(int64_t step, gmx::ArrayRef<const gmx::RVec> x, gmx::ForceWithVirial *force)
{

  //Get only the forces without virial
  rvec *f = as_rvec_array(force->force_.data());

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
    cvm::log(cvm::line_marker);
    cvm::log("colvarproxy_gromacs, step no. "+cvm::to_str(colvars->it)+"\n"+
             "Updating internal data.\n");
  }

  // backup applied forces if necessary to calculate total forces
  //if (total_force_requested)
  //  previous_atoms_new_colvar_forces = atoms_new_colvar_forces;

  // Zero the forces on the atoms, so that they can be accumulated by the colvars.
  for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++) {
    atoms_new_colvar_forces[i].x = atoms_new_colvar_forces[i].y = atoms_new_colvar_forces[i].z = 0.0;
  }

  // Get the atom positions from the Gromacs array.
  for (size_t i = 0; i < atoms_ids.size(); i++) {
    size_t aid = atoms_ids[i];
    if (aid >= gmx_atoms->nr) {
      cvm::fatal_error("Error: Atom index "+cvm::to_str(aid)+" not found in GROMACS data structure containing "+
                       cvm::to_str(gmx_atoms->nr)+" atoms");
    }
    atoms_positions[i] = cvm::rvector(x[aid][0], x[aid][1], x[aid][2]);
  }

  // Get total forces if required.
  if (total_force_requested && cvm::step_relative() > 0) {
     for (size_t i = 0; i < atoms_ids.size(); i++) {
       size_t aid = atoms_ids[i];
       // We already checked above that gmx_atoms->nr < aid.
       atoms_total_forces[i] = cvm::rvector(f[aid][0], f[aid][1], f[aid][2]);
     }
  }

  bias_energy = 0.0;
  // Call the collective variable module to fill atoms_new_colvar_forces
  if (colvars->calc() != COLVARS_OK) {
    cvm::fatal_error("");
  }

  // Pass the applied forces back to GROMACS.
  for (size_t i = 0; i < atoms_ids.size(); i++) {
    size_t aid = atoms_ids[i];
    // We already checked above that gmx_atoms->nr < aid.
    f[aid][0] += atoms_new_colvar_forces[i].x;
    f[aid][1] += atoms_new_colvar_forces[i].y;
    f[aid][2] += atoms_new_colvar_forces[i].z;
  }

  // We need to compute and update the virial like this (with virial as a 3x3 matrix):
  // matrix virial = compute_virial()
  // force->addVirialContribution(virial);

  return bias_energy;
}


// Pass restraint energy value for current timestep to MD engine
void colvarproxy_gromacs::add_energy (cvm::real energy)
{
  bias_energy += energy;
}

// **************** ATOMS ****************

int colvarproxy_gromacs::check_atom_id(int atom_number)
{
  // GROMACS uses zero-based arrays.
  int const aid = (atom_number-1);

  if (cvm::debug())
    log("Adding atom "+cvm::to_str(atom_number)+
        " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= gmx_atoms->nr) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  return aid;
}


int colvarproxy_gromacs::init_atom(int atom_number)
{
  // GROMACS uses zero-based arrays.
  int aid = atom_number-1;

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_ncopies[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);

  int const index = add_atom_slot(aid);
  update_atom_properties(index);
  return index;
}

void colvarproxy_gromacs::update_atom_properties(int index)
{

  // update mass
  double const mass = gmx_atoms->massT[atoms_ids[index]];
  if (mass <= 0.001) {
    this->log("Warning: near-zero mass for atom "+
              cvm::to_str(atoms_ids[index]+1)+
              "; expect unstable dynamics if you apply forces to it.\n");
  }
  atoms_masses[index] = mass;
  // update charge
  atoms_charges[index] = gmx_atoms->chargeA[atoms_ids[index]];
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
