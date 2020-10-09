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
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "colvarproxy_gromacs.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/domdec/ga2la.h"

#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/mdlib/broadcaststructs.h"

//************************************************************
// colvarproxy_gromacs
colvarproxy_gromacs::colvarproxy_gromacs() : colvarproxy() {}

// Colvars Initialization
void colvarproxy_gromacs::init(t_inputrec *ir, int64_t step,gmx_mtop_t *mtop,
                               const std::string &prefix,
                               gmx::ArrayRef<const std::string> filenames_config,
                               const std::string &filename_restart,
                               const t_commrec *cr,
                               const rvec x[]) {

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

  // GROMACS timestep
  timestep = ir->delta_t;
  // Retrieve the topology of all atoms
  gmx_atoms = gmx_mtop_global_atoms(mtop);

  // Read configuration file and set up the proxy only on the master node.
  if (MASTER(cr))
  {

    // initiate module: this object will be the communication proxy
    // colvarmodule pointer is only defined on the Master due to the static pointer to colvarproxy.
    colvars = new colvarmodule(this);

    version_int = get_version_from_string(COLVARPROXY_VERSION);

    if (cvm::debug()) {
      log("Initializing the colvars proxy object.\n");
    }

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
    }

    colvars->it = colvars->it_restart = step;

  } // end master


  // MPI initialisation

  // Initialise attributs for the MPI communication
  if(MASTER(cr)) {
    // Retrieve the number of colvar atoms
    nat = atoms_ids.size();
    // Copy their global indices
    ind = atoms_ids.data(); // This has to be updated if the vector is reallocated
  }


  if(PAR(cr)) {
    // Let the other nodes know the number of colvar atoms.
    block_bc(cr, nat);

    // Initialise atoms_new_colvar_forces on non-master nodes
    if(!MASTER(cr)) {
      atoms_new_colvar_forces.reserve(nat);
    }
  }

  snew(xa,         nat);
  snew(xa_ind,     nat);
  snew(xa_shifts,  nat);
  snew(xa_eshifts, nat);
  snew(xa_old,     nat);
  snew(f,          nat);

  // Prepare data

  // Save the original (whole) set of positions such that later the
  // molecule can always be made whole again
  if (MASTER(cr))
  {
    for (int i = 0; i < nat; i++)
    {
        int ii = ind[i];
        copy_rvec(x[ii], xa_old[i]);
    }
  }

  // Communicate initial coordinates and global indices to all processes
  if (PAR(cr))
  {
    nblock_bc(cr, nat, xa_old);
    snew_bc(cr, ind, nat);
    nblock_bc(cr, nat, ind);
  }

  // Serial Run
  if (!PAR(cr))
  {
    nat_loc = nat;
    nalloc_loc = nat;
    ind_loc = ind;

    // xa_ind[i] needs to be set to i for serial runs
    for (int i = 0; i < nat; i++)
    {
        xa_ind[i] = i;
    }
  }

  if (MASTER(cr) && cvm::debug()) {
    cvm::log ("atoms_ids = "+cvm::to_str (atoms_ids)+"\n");
    cvm::log ("atoms_ncopies = "+cvm::to_str (atoms_ncopies)+"\n");
    cvm::log ("positions = "+cvm::to_str (atoms_positions)+"\n");
    cvm::log ("total_forces = "+cvm::to_str (atoms_total_forces)+"\n");
    cvm::log ("atoms_new_colvar_forces = "+cvm::to_str (atoms_new_colvar_forces)+"\n");
    cvm::log (cvm::line_marker);
    log("done initializing the colvars proxy object.\n");
  }


} // End colvars initialization.


colvarproxy_gromacs::~colvarproxy_gromacs()
{
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
}

void colvarproxy_gromacs::finish(const t_commrec *cr)
{
  if(MASTER(cr)) {
    colvars->write_restart_file(output_prefix_str+".colvars.state");
    colvars->write_output_files();
  }
}

void colvarproxy_gromacs::set_temper(double temper)
{
  thermostat_temperature = temper;
}

// GROMACS uses nanometers and kJ/mol internally
cvm::real colvarproxy_gromacs::backend_angstrom_value() { return 0.1; }

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

void colvarproxy_gromacs::exit (std::string const gmx_unused &message)
{
  gmx_fatal(FARGS,"SUCCESS: %s\n", message.c_str());
}

int colvarproxy_gromacs::load_atoms (char const gmx_unused *filename, std::vector<cvm::atom> gmx_unused &atoms,
                                     std::string const gmx_unused &pdb_field, double const gmx_unused pdb_field_value)
{
  cvm::error("Selecting collective variable atoms "
		   "from a PDB file is currently not supported.\n");
  return COLVARS_NOT_IMPLEMENTED;
}

int colvarproxy_gromacs::load_coords (char const gmx_unused *filename, std::vector<cvm::atom_pos> gmx_unused &pos,
                                      const std::vector<int> gmx_unused &indices, std::string const gmx_unused &pdb_field_str,
                                      double const gmx_unused pdb_field_value)
{
  cvm::error("Selecting collective variable atoms "
		   "from a PDB file is currently not supported.\n");
  return COLVARS_NOT_IMPLEMENTED;
}

int colvarproxy_gromacs::set_unit_system(std::string const &units_in, bool /*colvars_defined*/)
{
  if (units_in != "gromacs") {
    cvm::error("Specified unit system \"" + units_in + "\" is unsupported in Gromacs. Supported units are \"gromacs\" (nm, kJ/mol).\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
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


void colvarproxy_gromacs::update_data(const t_commrec *cr, int64_t const step, t_pbc const &pbc, const matrix box, bool bNS)
{

  if (MASTER(cr)) {

    if(cvm::debug()) {
      cvm::log(cvm::line_marker);
      cvm::log("colvarproxy_gromacs, step no. "+cvm::to_str(colvars->it)+"\n"+
              "Updating internal data.\n");
    }

    // step update on master only due to the call of colvars pointer.
    if (first_timestep) {
      first_timestep = false;
    } else {
      // Use the time step number inherited from GROMACS
      if ( step - previous_gmx_step == 1 )
        colvars->it++;
      // Other cases?
    }
  } // end master

  gmx_pbc = pbc;
  gmx_box = box;
  gmx_bNS = bNS;

  previous_gmx_step = step;

  // Prepare data for MPI communication
  if(PAR(cr) && bNS) {
    dd_make_local_group_indices(cr->dd->ga2la, nat, ind, &nat_loc, &ind_loc, &nalloc_loc, xa_ind);
  }
}


void colvarproxy_gromacs::calculateForces(
                    const gmx::ForceProviderInput &forceProviderInput,
                    gmx::ForceProviderOutput      *forceProviderOutput)
{

  const t_commrec *cr           = &(forceProviderInput.cr_);
  const gmx::ArrayRef<const gmx::RVec> x  = forceProviderInput.x_;
  // Some gymnastics to coerce new data structures into old types
  const rvec *x_pointer          = &(x.data()->as_vec());


  // Eventually there needs to be an interface to update local data upon neighbor search
  // We could check if by chance all atoms are in one node, and skip communication
  communicate_group_positions(cr, xa, xa_shifts, xa_eshifts,
                              gmx_bNS, x_pointer, nat, nat_loc,
                              ind_loc, xa_ind, xa_old, gmx_box);

  // Communicate_group_positions takes care of removing shifts (unwrapping)
  // in single node jobs, communicate_group_positions() is efficient and adds no overhead

  if (MASTER(cr))
  {
    // On non-master nodes, jump directly to applying the forces

    // backup applied forces if necessary to calculate total forces
    //if (total_force_requested)
    //  previous_atoms_new_colvar_forces = atoms_new_colvar_forces;

    // Zero the forces on the atoms, so that they can be accumulated by the colvars.
    for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++) {
      atoms_new_colvar_forces[i].x = atoms_new_colvar_forces[i].y = atoms_new_colvar_forces[i].z = 0.0;
    }

    // Get the atom positions from the Gromacs array.
    for (size_t i = 0; i < atoms_ids.size(); i++) {
      atoms_positions[i] = cvm::rvector(xa[i][0], xa[i][1], xa[i][2]);
    }

    // // Get total forces if required.
    // if (total_force_requested && cvm::step_relative() > 0) {
    //   for (size_t i = 0; i < atoms_ids.size(); i++) {
    //     size_t aid = atoms_ids[i];
    //     // We already checked above that gmx_atoms->nr < aid.
    //     atoms_total_forces[i] = cvm::rvector(f[aid][0], f[aid][1], f[aid][2]);
    //   }
    // }

    bias_energy = 0.0;
    // Call the collective variable module to fill atoms_new_colvar_forces
    if (colvars->calc() != COLVARS_OK) {
      cvm::fatal_error("Error calling colvars->calc()\n");
    }

    // Copy the forces to a simpler array for broadcasting
    for (int i = 0; i < nat; i++)
    {
      f[i][0] = atoms_new_colvar_forces[i].x;
      f[i][1] = atoms_new_colvar_forces[i].y;
      f[i][2] = atoms_new_colvar_forces[i].z;
    }

    forceProviderOutput->enerd_.term[F_COM_PULL] += bias_energy;
  } // master node

  //Broadcast the forces to all the nodes
  if (PAR(cr))
  {
    nblock_bc(cr, nat, f);
  }


  const gmx::ArrayRef<gmx::RVec> &f_colvars = forceProviderOutput->forceWithVirial_.force_;

  // We need to compute and update the virial like this (with virial as a 3x3 matrix):
  // matrix virial = compute_virial()
  // force->addVirialContribution(virial);
  // virial is purely local, should be calculated where the atoms live

  // Pass the applied forces back to GROMACS
  // Parallel version
  for (int i = 0; i < nat; i++)
  {
      // j is the index in the "System group".
      int j = ind[i];

      // check if this is a local atom and find out locndx
      if (PAR(cr)) {
        const int *locndx = cr->dd->ga2la->findHome(j);
        if (locndx) {
          f_colvars[*locndx] += f[i];
        }
        // Do nothing if atom is not local
      } else { // Non MPI-parallel
        f_colvars[j] += f[i];
      }
  }

  return;
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

  if ( (aid < 0) || (aid >= gmx_atoms.nr) ) {
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

  if(aid < 0) {
    return INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);
  update_atom_properties(index);
  return index;
}

void colvarproxy_gromacs::update_atom_properties(int index)
{

  // update mass
  double const mass = gmx_atoms.atom[atoms_ids[index]].m;
  if (mass <= 0.001) {
    this->log("Warning: near-zero mass for atom "+
              cvm::to_str(atoms_ids[index]+1)+
              "; expect unstable dynamics if you apply forces to it.\n");
  }
  atoms_masses[index] = mass;
  // update charge
  atoms_charges[index] = gmx_atoms.atom[atoms_ids[index]].q;
}
