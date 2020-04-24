// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <cerrno>

#include <sstream>
#include <cstring>
#include <cstdio>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"
#include "colvaratoms.h"



colvarproxy_system::colvarproxy_system()
{
  angstrom_value = 0.0;
  total_force_requested = false;
  reset_pbc_lattice();
}


colvarproxy_system::~colvarproxy_system() {}


void colvarproxy_system::add_energy(cvm::real /* energy */) {}


void colvarproxy_system::request_total_force(bool yesno)
{
  if (yesno == true)
    cvm::error("Error: total forces are currently not implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
}


bool colvarproxy_system::total_forces_enabled() const
{
  return false;
}


bool colvarproxy_system::total_forces_same_step() const
{
  return false;
}


inline int round_to_integer(cvm::real x)
{
  return cvm::floor(x+0.5);
}


void colvarproxy_system::update_pbc_lattice()
{
  // Periodicity is assumed in all directions

  if (boundaries_type == boundaries_unsupported ||
      boundaries_type == boundaries_non_periodic) {
    cvm::error("Error: setting PBC lattice with unsupported boundaries.\n",
               BUG_ERROR);
    return;
  }

  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_y, unit_cell_z);
    reciprocal_cell_x = v/(v*unit_cell_x);
  }
  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_z, unit_cell_x);
    reciprocal_cell_y = v/(v*unit_cell_y);
  }
  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_x, unit_cell_y);
    reciprocal_cell_z = v/(v*unit_cell_z);
  }
}


void colvarproxy_system::reset_pbc_lattice()
{
  unit_cell_x.reset();
  unit_cell_y.reset();
  unit_cell_z.reset();
  reciprocal_cell_x.reset();
  reciprocal_cell_y.reset();
  reciprocal_cell_z.reset();
}


cvm::rvector colvarproxy_system::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2)
  const
{
  if (boundaries_type == boundaries_unsupported) {
    cvm::error("Error: unsupported boundary conditions.\n", INPUT_ERROR);
  }

  cvm::rvector diff = (pos2 - pos1);

  if (boundaries_type == boundaries_non_periodic) return diff;

  cvm::real const x_shift = round_to_integer(reciprocal_cell_x*diff);
  cvm::real const y_shift = round_to_integer(reciprocal_cell_y*diff);
  cvm::real const z_shift = round_to_integer(reciprocal_cell_z*diff);

  diff.x -= x_shift*unit_cell_x.x + y_shift*unit_cell_y.x +
    z_shift*unit_cell_z.x;
  diff.y -= x_shift*unit_cell_x.y + y_shift*unit_cell_y.y +
    z_shift*unit_cell_z.y;
  diff.z -= x_shift*unit_cell_x.z + y_shift*unit_cell_y.z +
    z_shift*unit_cell_z.z;

  return diff;
}


int colvarproxy_system::get_molid(int &)
{
  cvm::error("Error: only VMD allows the use of multiple \"molecules\", "
             "i.e. multiple molecular systems.", COLVARS_NOT_IMPLEMENTED);
  return -1;
}



colvarproxy_atoms::colvarproxy_atoms()
{
  updated_masses_ = updated_charges_ = false;
}


colvarproxy_atoms::~colvarproxy_atoms()
{
  reset();
}


int colvarproxy_atoms::reset()
{
  atoms_ids.clear();
  atoms_ncopies.clear();
  atoms_masses.clear();
  atoms_charges.clear();
  atoms_positions.clear();
  atoms_total_forces.clear();
  atoms_new_colvar_forces.clear();
  return COLVARS_OK;
}


int colvarproxy_atoms::add_atom_slot(int atom_id)
{
  atoms_ids.push_back(atom_id);
  atoms_ncopies.push_back(1);
  atoms_masses.push_back(1.0);
  atoms_charges.push_back(0.0);
  atoms_positions.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atoms_total_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atoms_new_colvar_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  return (atoms_ids.size() - 1);
}


int colvarproxy_atoms::init_atom(cvm::residue_id const & /* residue */,
                                 std::string const     & /* atom_name */,
                                 std::string const     & /* segment_id */)
{
  cvm::error("Error: initializing an atom by name and residue number is currently not supported.\n",
             COLVARS_NOT_IMPLEMENTED);
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_atoms::check_atom_id(cvm::residue_id const &residue,
                                     std::string const     &atom_name,
                                     std::string const     &segment_id)
{
  colvarproxy_atoms::init_atom(residue, atom_name, segment_id);
  return COLVARS_NOT_IMPLEMENTED;
}


void colvarproxy_atoms::clear_atom(int index)
{
  if (((size_t) index) >= atoms_ids.size()) {
    cvm::error("Error: trying to disable an atom that was not previously requested.\n",
               INPUT_ERROR);
  }
  if (atoms_ncopies[index] > 0) {
    atoms_ncopies[index] -= 1;
  }
}


int colvarproxy_atoms::load_atoms(char const * /* filename */,
                                  cvm::atom_group & /* atoms */,
                                  std::string const & /* pdb_field */,
                                  double)
{
  return cvm::error("Error: loading atom identifiers from a file "
                    "is currently not implemented.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_atoms::load_coords(char const * /* filename */,
                                   std::vector<cvm::atom_pos> & /* pos */,
                                   std::vector<int> const & /* sorted_ids */,
                                   std::string const & /* pdb_field */,
                                   double)
{
  return cvm::error("Error: loading atomic coordinates from a file "
                    "is currently not implemented.\n",
                    COLVARS_NOT_IMPLEMENTED);
}



colvarproxy_atom_groups::colvarproxy_atom_groups() {}


colvarproxy_atom_groups::~colvarproxy_atom_groups()
{
  reset();
}


int colvarproxy_atom_groups::reset()
{
  atom_groups_ids.clear();
  atom_groups_ncopies.clear();
  atom_groups_masses.clear();
  atom_groups_charges.clear();
  atom_groups_coms.clear();
  atom_groups_total_forces.clear();
  atom_groups_new_colvar_forces.clear();
  return COLVARS_OK;
}


int colvarproxy_atom_groups::add_atom_group_slot(int atom_group_id)
{
  atom_groups_ids.push_back(atom_group_id);
  atom_groups_ncopies.push_back(1);
  atom_groups_masses.push_back(1.0);
  atom_groups_charges.push_back(0.0);
  atom_groups_coms.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atom_groups_total_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atom_groups_new_colvar_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  return (atom_groups_ids.size() - 1);
}


int colvarproxy_atom_groups::scalable_group_coms()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_atom_groups::init_atom_group(std::vector<int> const & /* atoms_ids */)
{
  cvm::error("Error: initializing a group outside of the Colvars module "
             "is currently not supported.\n",
             COLVARS_NOT_IMPLEMENTED);
  return COLVARS_NOT_IMPLEMENTED;
}


void colvarproxy_atom_groups::clear_atom_group(int index)
{
  if (((size_t) index) >= atom_groups_ids.size()) {
    cvm::error("Error: trying to disable an atom group "
               "that was not previously requested.\n",
               INPUT_ERROR);
  }
  if (atom_groups_ncopies[index] > 0) {
    atom_groups_ncopies[index] -= 1;
  }
}



colvarproxy_smp::colvarproxy_smp()
{
  b_smp_active = true; // May be disabled by user option
  omp_lock_state = NULL;
#if defined(_OPENMP)
  if (smp_thread_id() == 0) {
    omp_lock_state = reinterpret_cast<void *>(new omp_lock_t);
    omp_init_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state));
  }
#endif
}


colvarproxy_smp::~colvarproxy_smp()
{
#if defined(_OPENMP)
  if (smp_thread_id() == 0) {
    if (omp_lock_state) {
      delete reinterpret_cast<omp_lock_t *>(omp_lock_state);
    }
  }
#endif
}


int colvarproxy_smp::smp_enabled()
{
#if defined(_OPENMP)
  if (b_smp_active) {
    return COLVARS_OK;
  }
  return COLVARS_ERROR;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_colvars_loop()
{
#if defined(_OPENMP)
  colvarmodule *cv = cvm::main();
  colvarproxy *proxy = cv->proxy;
#pragma omp parallel for
  for (size_t i = 0; i < cv->variables_active_smp()->size(); i++) {
    colvar *x = (*(cv->variables_active_smp()))[i];
    int x_item = (*(cv->variables_active_smp_items()))[i];
    if (cvm::debug()) {
      cvm::log("["+cvm::to_str(proxy->smp_thread_id())+"/"+
               cvm::to_str(proxy->smp_num_threads())+
               "]: calc_colvars_items_smp(), i = "+cvm::to_str(i)+", cv = "+
               x->name+", cvc = "+cvm::to_str(x_item)+"\n");
    }
    x->calc_cvcs(x_item, 1);
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_biases_loop()
{
#if defined(_OPENMP)
  colvarmodule *cv = cvm::main();
#pragma omp parallel
  {
#pragma omp for
    for (size_t i = 0; i < cv->biases_active()->size(); i++) {
      colvarbias *b = (*(cv->biases_active()))[i];
      if (cvm::debug()) {
        cvm::log("Calculating bias \""+b->name+"\" on thread "+
                 cvm::to_str(smp_thread_id())+"\n");
      }
      b->update();
    }
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_biases_script_loop()
{
#if defined(_OPENMP)
  colvarmodule *cv = cvm::main();
#pragma omp parallel
  {
#pragma omp single nowait
    {
      cv->calc_scripted_forces();
    }
#pragma omp for
    for (size_t i = 0; i < cv->biases_active()->size(); i++) {
      colvarbias *b = (*(cv->biases_active()))[i];
      if (cvm::debug()) {
        cvm::log("Calculating bias \""+b->name+"\" on thread "+
                 cvm::to_str(smp_thread_id())+"\n");
      }
      b->update();
    }
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}




int colvarproxy_smp::smp_thread_id()
{
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_num_threads()
{
#if defined(_OPENMP)
  return omp_get_max_threads();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_lock()
{
#if defined(_OPENMP)
  omp_set_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state));
#endif
  return COLVARS_OK;
}


int colvarproxy_smp::smp_trylock()
{
#if defined(_OPENMP)
  return omp_test_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state)) ?
    COLVARS_OK : COLVARS_ERROR;
#else
  return COLVARS_OK;
#endif
}


int colvarproxy_smp::smp_unlock()
{
#if defined(_OPENMP)
  omp_unset_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state));
#endif
  return COLVARS_OK;
}



colvarproxy_script::colvarproxy_script()
{
  script = NULL;
}


colvarproxy_script::~colvarproxy_script() {}


char const *colvarproxy_script::script_obj_to_str(unsigned char *obj)
{
  cvm::error("Error: trying to print a script object without a scripting "
             "language interface.\n", BUG_ERROR);
  return reinterpret_cast<char *>(obj);
}


std::vector<std::string> colvarproxy_script::script_obj_to_str_vector(unsigned char * /* obj */)
{
  cvm::error("Error: trying to print a script object without a scripting "
             "language interface.\n", BUG_ERROR);
  return std::vector<std::string>();
}


int colvarproxy_script::run_force_callback()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_script::run_colvar_callback(std::string const & /* name */,
					    std::vector<const colvarvalue *> const & /* cvcs */,
					    colvarvalue & /* value */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_script::run_colvar_gradient_callback(std::string const & /* name */,
						     std::vector<const colvarvalue *> const & /* cvcs */,
						     std::vector<cvm::matrix2d<cvm::real> > & /* gradient */)
{
  return COLVARS_NOT_IMPLEMENTED;
}



colvarproxy_io::colvarproxy_io()
{
  input_buffer_ = NULL;
}


colvarproxy_io::~colvarproxy_io() {}


int colvarproxy_io::get_frame(long int&)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::set_frame(long int)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::backup_file(char const * /* filename */)
{
  // TODO implement this using rename_file()
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::remove_file(char const *filename)
{
  int error_code = COLVARS_OK;
#if defined(WIN32) && !defined(__CYGWIN__)
  // Because the file may be open by other processes, rename it to filename.old
  std::string const renamed_file(std::string(filename)+".old");
  // It may still be there from an interrupted run, so remove it to be safe
  std::remove(renamed_file.c_str());
  int rename_exit_code = 0;
  while ((rename_exit_code = std::rename(filename,
                                         renamed_file.c_str())) != 0) {
    if (errno == EINTR) continue;
    error_code |= FILE_ERROR;
    break;
  }
  // Ask to remove filename.old, but ignore any errors raised
  std::remove(renamed_file.c_str());
#else
  if (std::remove(filename)) {
    if (errno != ENOENT) {
      error_code |= FILE_ERROR;
    }
  }
#endif
  if (error_code != COLVARS_OK) {
    return cvm::error("Error: in removing file \""+std::string(filename)+
                      "\".\n.",
                      error_code);
  }
  return COLVARS_OK;
}


int colvarproxy_io::rename_file(char const *filename, char const *newfilename)
{
  int error_code = COLVARS_OK;
#if defined(WIN32) && !defined(__CYGWIN__)
  // On straight Windows, must remove the destination before renaming it
  error_code |= remove_file(newfilename);
#endif
  int rename_exit_code = 0;
  while ((rename_exit_code = std::rename(filename, newfilename)) != 0) {
    if (errno == EINTR) continue;
    // Call log() instead of error to allow the next try
    cvm::log("Error: in renaming file \""+std::string(filename)+"\" to \""+
             std::string(newfilename)+"\".\n.");
    error_code |= FILE_ERROR;
    if (errno == EXDEV) continue;
    break;
  }
  return rename_exit_code ? error_code : COLVARS_OK;
}



colvarproxy::colvarproxy()
{
  colvars = NULL;
  b_simulation_running = true;
  b_delete_requested = false;
}


colvarproxy::~colvarproxy()
{
  close_files();
}


int colvarproxy::close_files()
{
  if (smp_enabled() == COLVARS_OK && smp_thread_id() > 0) {
    // Nothing to do on non-master threads
    return COLVARS_OK;
  }
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    ((std::ofstream *) (*osi))->close();
    delete *osi;
  }
  output_files.clear();
  output_stream_names.clear();
  return COLVARS_OK;
}


int colvarproxy::reset()
{
  int error_code = COLVARS_OK;
  error_code |= colvarproxy_atoms::reset();
  error_code |= colvarproxy_atom_groups::reset();
  return error_code;
}


int colvarproxy::request_deletion()
{
  return cvm::error("Error: \"delete\" command is only available in VMD; "
                    "please use \"reset\" instead.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy::setup()
{
  return COLVARS_OK;
}


int colvarproxy::update_input()
{
  return COLVARS_OK;
}


int colvarproxy::update_output()
{
  return COLVARS_OK;
}


int colvarproxy::post_run()
{
  int error_code = COLVARS_OK;
  if (colvars->output_prefix().size()) {
    error_code |= colvars->write_restart_file(cvm::output_prefix()+".colvars.state");
    error_code |= colvars->write_output_files();
  }
  error_code |= flush_output_streams();
  return error_code;
}


void colvarproxy::add_error_msg(std::string const &message)
{
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line)) {
    error_output += line+"\n";
  }
}


void colvarproxy::clear_error_msgs()
{
  error_output.clear();
}


std::string const & colvarproxy::get_error_msgs()
{
  return error_output;
}


int colvarproxy::get_version_from_string(char const *version_string)
{
  std::string const v(version_string);
  std::istringstream is(v.substr(0, 4) + v.substr(5, 2) + v.substr(8, 2));
  int newint;
  is >> newint;
  return newint;
}


void colvarproxy::smp_stream_error()
{
  cvm::error("Error: trying to access an output stream from a "
             "multi-threaded region (bug).  For a quick workaround, use "
             "\"smp off\" in the Colvars config.\n", BUG_ERROR);
}


std::ostream * colvarproxy::output_stream(std::string const &output_name,
                                          std::ios_base::openmode mode)
{
  if (cvm::debug()) {
    cvm::log("Using colvarproxy::output_stream()\n");
  }

  std::ostream *os = get_output_stream(output_name);
  if (os != NULL) return os;

  if (!(mode & (std::ios_base::app | std::ios_base::ate))) {
    backup_file(output_name);
  }
  std::ofstream *osf = new std::ofstream(output_name.c_str(), mode);
  if (!osf->is_open()) {
    cvm::error("Error: cannot write to file/channel \""+output_name+"\".\n",
               FILE_ERROR);
    return NULL;
  }
  output_stream_names.push_back(output_name);
  output_files.push_back(osf);
  return osf;
}


std::ostream *colvarproxy::get_output_stream(std::string const &output_name)
{
  if (smp_enabled() == COLVARS_OK) {
    if (smp_thread_id() > 0) smp_stream_error();
  }
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osni == output_name) {
      return *osi;
    }
  }
  return NULL;
}



int colvarproxy::flush_output_stream(std::ostream *os)
{
  if (smp_enabled() == COLVARS_OK) {
    if (smp_thread_id() > 0) smp_stream_error();
  }
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osi == os) {
      ((std::ofstream *) (*osi))->flush();
      return COLVARS_OK;
    }
  }
  return cvm::error("Error: trying to flush an output file/channel "
                    "that wasn't open.\n", BUG_ERROR);
}


int colvarproxy::flush_output_streams()
{
  if (smp_enabled() == COLVARS_OK && smp_thread_id() > 0)
    return COLVARS_OK;

  std::list<std::ostream *>::iterator osi  = output_files.begin();
  for ( ; osi != output_files.end(); osi++) {
    ((std::ofstream *) (*osi))->flush();
  }
  return COLVARS_OK;
}


int colvarproxy::close_output_stream(std::string const &output_name)
{
  if (smp_enabled() == COLVARS_OK) {
    if (smp_thread_id() > 0) smp_stream_error();
  }
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osni == output_name) {
      ((std::ofstream *) (*osi))->close();
      delete *osi;
      output_files.erase(osi);
      output_stream_names.erase(osni);
      return COLVARS_OK;
    }
  }
  return cvm::error("Error: trying to close an output file/channel "
                    "that wasn't open.\n", BUG_ERROR);
}
