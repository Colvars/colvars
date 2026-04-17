// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <fstream>
#include <list>
#include <utility>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvarscript.h"
#include "colvarmodule_utils.h"

#include <algorithm>

#ifdef COLVARS_TCL
#include <tcl.h>
#endif

int colvarproxy::reset_atoms()
{
  atoms_ids.clear();
  atoms_refcount.clear();
  atoms_masses.clear();
  atoms_charges.clear();
  atoms_positions.clear();
  atoms_total_forces.clear();
  atoms_new_colvar_forces.clear();
  return COLVARS_OK;
}


int colvarproxy::add_atom_slot(int atom_id)
{
  atoms_ids.push_back(atom_id);
  atoms_refcount.push_back(1);
  atoms_masses.push_back(1.0);
  atoms_charges.push_back(0.0);
  atoms_positions.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atoms_total_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atoms_new_colvar_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  modified_atom_list_ = true;
  return (atoms_ids.size() - 1);
}


int colvarproxy::init_atom(int /* atom_number */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy::check_atom_id(int /* atom_number */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy::check_atom_name_selections_available()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy::init_atom(cvm::residue_id const & /* residue */,
                                 std::string const     & /* atom_name */,
                                 std::string const     & /* segment_id */)
{
  cvm::error_static(cvmodule, "Error: initializing an atom by name and residue number is currently not supported.\n",
             COLVARS_NOT_IMPLEMENTED);
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy::check_atom_id(cvm::residue_id const &residue,
                                     std::string const     &atom_name,
                                     std::string const     &segment_id)
{
  init_atom(residue, atom_name, segment_id);
  return COLVARS_NOT_IMPLEMENTED;
}


void colvarproxy::clear_atom(int index)
{
  if (((size_t) index) >= atoms_ids.size()) {
    cvm::error_static(cvmodule, "Error: trying to disable an atom that was not previously requested.\n",
               COLVARS_INPUT_ERROR);
  }
  if (atoms_refcount[index] > 0) {
    atoms_refcount[index] -= 1;
  }
}


size_t colvarproxy::get_num_active_atoms() const
{
  size_t result = 0;
  for (size_t i = 0; i < atoms_refcount.size(); i++) {
    if (atoms_refcount[i] > 0) result++;
  }
  return result;
}


void colvarproxy::compute_rms_atoms_applied_force()
{
  atoms_rms_applied_force_ =
    compute_norm2_stats<decltype(atoms_new_colvar_forces), 0, false>(atoms_new_colvar_forces);
}


void colvarproxy::compute_max_atoms_applied_force()
{
  int minmax_index = -1;
  size_t const n_atoms_ids = atoms_ids.size();
  if ((n_atoms_ids > 0) && (n_atoms_ids == atoms_new_colvar_forces.size())) {
    atoms_max_applied_force_ =
      compute_norm2_stats<decltype(atoms_new_colvar_forces), 1, true>(atoms_new_colvar_forces,
                                                 &minmax_index);
    if (minmax_index >= 0) {
      atoms_max_applied_force_id_ = atoms_ids[minmax_index];
    } else {
      atoms_max_applied_force_id_ = -1;
    }
  } else {
    atoms_max_applied_force_ =
      compute_norm2_stats<decltype(atoms_new_colvar_forces), 1, false>(atoms_new_colvar_forces);
    atoms_max_applied_force_id_ = -1;
  }
}


int colvarproxy::reset_atom_groups()
{
  atom_groups_ids.clear();
  atom_groups_refcount.clear();
  atom_groups_masses.clear();
  atom_groups_charges.clear();
  atom_groups_coms.clear();
  atom_groups_total_forces.clear();
  atom_groups_new_colvar_forces.clear();
  return COLVARS_OK;
}


int colvarproxy::add_atom_group_slot(int atom_group_id)
{
  atom_groups_ids.push_back(atom_group_id);
  atom_groups_refcount.push_back(1);
  atom_groups_masses.push_back(1.0);
  atom_groups_charges.push_back(0.0);
  atom_groups_coms.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atom_groups_total_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atom_groups_new_colvar_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  return (atom_groups_ids.size() - 1);
}


int colvarproxy::scalable_group_coms()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy::init_atom_group(std::vector<int> const & /* atoms_ids */)
{
  cvm::error_static(cvmodule, "Error: initializing a group outside of the Colvars module "
             "is currently not supported.\n",
             COLVARS_NOT_IMPLEMENTED);
  return COLVARS_NOT_IMPLEMENTED;
}


void colvarproxy::clear_atom_group(int index)
{
  if (((size_t) index) >= atom_groups_ids.size()) {
    cvm::error_static(cvmodule, "Error: trying to disable an atom group "
               "that was not previously requested.\n",
               COLVARS_INPUT_ERROR);
  }
  if (atom_groups_refcount[index] > 0) {
    atom_groups_refcount[index] -= 1;
  }
}


size_t colvarproxy::get_num_active_atom_groups() const
{
  size_t result = 0;
  for (size_t i = 0; i < atom_groups_refcount.size(); i++) {
    if (atom_groups_refcount[i] > 0) result++;
  }
  return result;
}


void colvarproxy::compute_rms_atom_groups_applied_force()
{
  atom_groups_rms_applied_force_ =
    compute_norm2_stats<decltype(atom_groups_new_colvar_forces), 0, false>(atom_groups_new_colvar_forces);
}


void colvarproxy::compute_max_atom_groups_applied_force()
{
  atom_groups_max_applied_force_ =
    compute_norm2_stats<decltype(atom_groups_new_colvar_forces), 1, false>(atom_groups_new_colvar_forces);
}



colvarproxy::smp_mode_t colvarproxy::get_smp_mode() const {
  return smp_mode;
}

std::vector<colvarproxy::smp_mode_t> colvarproxy::get_available_smp_modes() const {
  std::vector<colvarproxy::smp_mode_t> modes;
#if defined(_OPENMP)
  modes.push_back(colvarproxy::smp_mode_t::cvcs);
  modes.push_back(colvarproxy::smp_mode_t::inner_loop);
#endif
  modes.push_back(colvarproxy::smp_mode_t::none);
  return modes;
}

colvarproxy::smp_mode_t colvarproxy::get_preferred_smp_mode() const {
  return get_available_smp_modes()[0];
}

int colvarproxy::set_smp_mode(smp_mode_t mode) {
  std::vector<colvarproxy::smp_mode_t> available_modes = get_available_smp_modes();
  auto it = std::find(available_modes.begin(), available_modes.end(), mode);
  if (it != available_modes.end()) {
    smp_mode = *it;
    return COLVARS_OK;
  } else {
    return COLVARS_NOT_IMPLEMENTED;
  }
}


int colvarproxy::smp_loop(int n_items, std::function<int (int)> const &worker)
{
  int error_code = COLVARS_OK;
#if defined(_OPENMP)
  if (cvmodule) cvmodule->increase_depth();
#pragma omp parallel for
  for (int i = 0; i < n_items; i++) {
    int const retcode = worker(i);
#pragma omp atomic
    error_code |= retcode;
  }
  if (cvmodule) cvmodule->decrease_depth();
#else
  error_code |= COLVARS_NOT_IMPLEMENTED;
#endif
  return error_code;
}


int colvarproxy::smp_biases_loop()
{
  if (!cvmodule) return COLVARS_BUG_ERROR;
#if defined(_OPENMP)
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < static_cast<int>(cvmodule->biases_active()->size()); i++) {
      colvarbias *b = (*(cvmodule->biases_active()))[i];
      if (cvmodule->debug()) {
        cvmodule->log("Calculating bias \""+b->name+"\" on thread "+
                 cvm::to_str(smp_thread_id())+"\n");
      }
      b->update();
    }
  }
  return cvmodule->get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy::smp_biases_script_loop()
{
#if defined(_OPENMP)
  if (!cvmodule) return COLVARS_BUG_ERROR;
#pragma omp parallel
  {
#pragma omp single nowait
    {
      cvmodule->calc_scripted_forces();
    }
#pragma omp for
    for (int i = 0; i < static_cast<int>(cvmodule->biases_active()->size()); i++) {
      colvarbias *b = (*(cvmodule->biases_active()))[i];
      if (cvmodule->debug()) {
        cvmodule->log("Calculating bias \""+b->name+"\" on thread "+
                 cvm::to_str(smp_thread_id())+"\n");
      }
      b->update();
    }
  }
  return cvmodule->get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}

int colvarproxy::smp_thread_id()
{
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return -1;
#endif
}


int colvarproxy::smp_num_threads()
{
#if defined(_OPENMP)
  return omp_get_max_threads();
#else
  return -1;
#endif
}


int colvarproxy::smp_lock()
{
#if defined(_OPENMP)
  omp_set_lock(omp_lock_state);
#endif
  return COLVARS_OK;
}


int colvarproxy::smp_trylock()
{
#if defined(_OPENMP)
  return omp_test_lock(omp_lock_state) ? COLVARS_OK : COLVARS_ERROR;
#else
  return COLVARS_OK;
#endif
}


int colvarproxy::smp_unlock()
{
#if defined(_OPENMP)
  omp_unset_lock(omp_lock_state);
#endif
  return COLVARS_OK;
}




int colvarproxy::run_force_callback()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy::run_colvar_callback(std::string const & /* name */,
                                            std::vector<const colvarvalue *> const & /* cvcs */,
                                            colvarvalue & /* value */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy::run_colvar_gradient_callback(std::string const & /* name */,
                                                     std::vector<const colvarvalue *> const & /* cvcs */,
                                                     std::vector<cvm::matrix2d<cvm::real> > & /* gradient */)
{
  return COLVARS_NOT_IMPLEMENTED;
}



colvarproxy::colvarproxy()
{
  cvmodule = NULL;
  // colvarproxy_system
  angstrom_value_ = 0.0;
  kcal_mol_value_ = 0.0;
  timestep_ = 1.0;
  target_temperature_ = 0.0;
  boltzmann_ = 0.001987191; // Default: kcal/mol/K
  boundaries_type = boundaries_unsupported;
  total_force_requested = false;
  indirect_lambda_biasing_force = 0.0;
  cached_alch_lambda_changed = false;
  cached_alch_lambda = -1.0;
  reset_pbc_lattice();
  // colvarproxy_atoms
  atoms_rms_applied_force_ = atoms_max_applied_force_ = 0.0;
  atoms_max_applied_force_id_ = -1;
  modified_atom_list_ = false;
  updated_masses_ = updated_charges_ = false;
  // colvarproxy_atom_groups
  atom_groups_rms_applied_force_ = atom_groups_max_applied_force_ = 0.0;
  // colvarproxy_volmaps
  volmaps_rms_applied_force_ = volmaps_max_applied_force_ = 0.0;
  // colvarproxy_smp
  smp_mode = smp_mode_t::none; // May be disabled by user option
  omp_lock_state = NULL;
#if defined(_OPENMP)
  smp_mode = smp_mode_t::cvcs;
  if (omp_get_thread_num() == 0) {
    omp_lock_state = new omp_lock_t;
    omp_init_lock(omp_lock_state);
  }
#endif
  // colvarproxy_replicas
#ifdef COLVARS_MPI
  replicas_mpi_comm = MPI_COMM_NULL;
#endif
  // colvarproxy_tcl
  tcl_interp_ = nullptr;
  // colvarproxy_io
  restart_frequency_engine = 0;
  input_stream_error_ = new std::istringstream();
  input_stream_error_->setstate(std::ios::badbit);
  output_stream_error_ = new std::ostringstream();
  output_stream_error_->setstate(std::ios::badbit);
  // colvarproxy_gpu
  support_gpu = false;
  // By default, simulation engines allow to immediately request atoms
  engine_ready_ = true;
  b_simulation_running = true;
  b_simulation_continuing = false;
  b_delete_requested = false;
  version_int = -1;
  features_hash = 0;
  config_queue_ = reinterpret_cast<void *>(new std::list<std::pair<std::string, std::string> >);
}


colvarproxy::~colvarproxy()
{
  // Destruct colvarproxy_script
  if (script) {
    delete script;
    script = nullptr;
  }
  // The following dtor will be called by cvm::reset from the dtor of colvarmodule??
  // // Destruct colvarproxy_volmaps
  // reset_volmaps();
  // // Destruct colvarproxy_atom_groups
  // reset_atom_groups();
  // // Destruct colvarproxy_atoms
  // reset_atoms();
  if (cvmodule) {
    delete cvmodule;
    cvmodule = nullptr;
  }
  // TODO: I guess that "delete cvmodule" will need to flush the streams, so close them after destructing cvmodule.
  // Destruct colvarproxy_io
  close_input_streams();
  close_output_streams();
  if (input_stream_error_) {
    delete input_stream_error_;
  }
  if (output_stream_error_) {
    delete output_stream_error_;
  }
  delete reinterpret_cast<std::list<std::pair<std::string, std::string> > *>(config_queue_);
#if defined(_OPENMP)
  if (omp_get_thread_num() == 0) {
    if (omp_lock_state) {
      delete omp_lock_state;
    }
  }
#endif
}


bool colvarproxy::io_available()
{
  return ((get_smp_mode() != smp_mode_t::none) && smp_thread_id() == 0) ||
         (get_smp_mode() == smp_mode_t::none) ||
         (get_smp_mode() == smp_mode_t::gpu);
}


int colvarproxy::reset()
{
  if (cvm::debug()) {
    cvmodule->log("colvarproxy::reset()\n");
  }
  int error_code = COLVARS_OK;
  error_code |= reset_atoms();
  error_code |= reset_atom_groups();
  error_code |= reset_volmaps();
  total_force_requested = false;
  return error_code;
}


int colvarproxy::request_deletion()
{
  return cvmodule->error("Error: \"delete\" command is only available in VMD; "
                    "please use \"reset\" instead.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


void colvarproxy::add_config(std::string const &cmd, std::string const &conf)
{
  reinterpret_cast<std::list<std::pair<std::string, std::string> > *>(config_queue_)->push_back(std::make_pair(cmd, conf));
}


int colvarproxy::setup()
{
  return COLVARS_OK;
}


int colvarproxy::parse_module_config()
{
  int error_code = COLVARS_OK;
  // Read any configuration queued up for Colvars
  std::list<std::pair<std::string, std::string> > *config_queue = reinterpret_cast<std::list<std::pair<std::string, std::string> > *>(config_queue_);
  while (config_queue->size() > 0) {
    std::pair<std::string, std::string> const &p = config_queue->front();
    if (p.first == "config") {
      error_code |= cvmodule->read_config_string(p.second);
    } else if (p.first == "configfile") {
      error_code |= cvmodule->read_config_file(p.second.c_str());
    } else {
      error_code |= cvmodule->error(std::string("Error: invalid keyword \"") +
                               p.first +
                               std::string("\" in colvarproxy::setup()\n"),
                               COLVARS_BUG_ERROR);
    }
    config_queue->pop_front();
  }
  return error_code;
}

int colvarproxy::load_atoms_pdb(char const * /* filename */,
                                cvm::atom_group & /* atoms */,
                                std::string const & /* pdb_field */,
                                double /* pdb_field_value */)
{
  return cvmodule->error(
      "Error: loading atom indices from a PDB file is currently not implemented in " +
          engine_name() + ".\n",
      COLVARS_NOT_IMPLEMENTED);
}

int colvarproxy::load_coords_pdb(char const * /* filename */,
                                 std::vector<cvm::atom_pos> & /* pos */,
                                 std::vector<int> const & /* sorted_ids */,
                                 std::string const & /* pdb_field */,
                                 double /* pdb_field_value */)
{
  return cvmodule->error(
      "Error: loading atomic coordinates from a PDB file is currently not implemented in " +
          engine_name() + ".\n",
      COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy::update_input()
{
  return COLVARS_OK;
}


int colvarproxy::update_output()
{
  return COLVARS_OK;
}


int colvarproxy::end_of_step()
{
  // Disable flags that Colvars doesn't need any more
  updated_masses_ = updated_charges_ = false;

  // Compute force statistics
  compute_rms_atoms_applied_force();
  compute_max_atoms_applied_force();
  compute_rms_atom_groups_applied_force();
  compute_max_atom_groups_applied_force();
  compute_rms_volmaps_applied_force();
  compute_max_volmaps_applied_force();

  if (cached_alch_lambda_changed) {
    send_alch_lambda();
    cached_alch_lambda_changed = false;
  }
  return COLVARS_OK;
}


int colvarproxy::post_run()
{
  int error_code = COLVARS_OK;
  if (cvmodule->output_prefix().size()) {
    error_code |= cvmodule->write_restart_file(cvmodule->output_prefix()+".colvars.state");
    error_code |= cvmodule->write_output_files();
  }
  error_code |= flush_output_streams();
  return error_code;
}


void colvarproxy::set_total_forces_invalid()
{
  std::fill(atoms_total_forces.begin(), atoms_total_forces.end(), cvm::rvector(0.0, 0.0, 0.0));
  std::fill(atom_groups_total_forces.begin(), atom_groups_total_forces.end(),
            cvm::rvector(0.0, 0.0, 0.0));
  total_forces_valid_ = false;
}


void colvarproxy::set_total_forces_valid()
{
  total_forces_valid_ = true;
}


void colvarproxy::print_input_atomic_data()
{
  cvmodule->log(cvm::line_marker);

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atoms_ids[size = "+cvm::to_str(atoms_ids.size())+
           "] = "+cvm::to_str(atoms_ids)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atoms_refcount[size = "+cvm::to_str(atoms_refcount.size())+
           "] = "+cvm::to_str(atoms_refcount)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atoms_masses[size = "+cvm::to_str(atoms_masses.size())+
           "] = "+cvm::to_str(atoms_masses)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atoms_charges[size = "+cvm::to_str(atoms_charges.size())+
           "] = "+cvm::to_str(atoms_charges)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atoms_positions[size = "+cvm::to_str(atoms_positions.size())+
           "] = "+cvm::to_str(atoms_positions,
                              cvmodule->cv_width,
                              cvmodule->cv_prec)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atoms_total_forces[size = "+
           cvm::to_str(atoms_total_forces.size())+
           "] = "+cvm::to_str(atoms_total_forces,
                              cvmodule->cv_width,
                              cvmodule->cv_prec)+"\n");

  cvmodule->log(cvm::line_marker);

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atom_groups_ids[size = "+cvm::to_str(atom_groups_ids.size())+
           "] = "+cvm::to_str(atom_groups_ids)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atom_groups_refcount[size = "+
           cvm::to_str(atom_groups_refcount.size())+
           "] = "+cvm::to_str(atom_groups_refcount)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atom_groups_masses[size = "+
           cvm::to_str(atom_groups_masses.size())+
           "] = "+cvm::to_str(atom_groups_masses)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atom_groups_charges[size = "+
           cvm::to_str(atom_groups_charges.size())+
           "] = "+cvm::to_str(atom_groups_charges)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atom_groups_coms[size = "+
           cvm::to_str(atom_groups_coms.size())+
           "] = "+cvm::to_str(atom_groups_coms,
                              cvmodule->cv_width,
                              cvmodule->cv_prec)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atom_groups_total_forces[size = "+
           cvm::to_str(atom_groups_total_forces.size())+
           "] = "+cvm::to_str(atom_groups_total_forces,
                              cvmodule->cv_width,
                              cvmodule->cv_prec)+"\n");

  cvmodule->log(cvm::line_marker);

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "volmaps_ids[size = "+cvm::to_str(volmaps_ids.size())+
           "] = "+cvm::to_str(volmaps_ids)+"\n");

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "volmaps_values[size = "+cvm::to_str(volmaps_values.size())+
           "] = "+cvm::to_str(volmaps_values)+"\n");

  cvmodule->log(cvm::line_marker);
}


void colvarproxy::print_output_atomic_data()
{
  cvmodule->log(cvm::line_marker);
  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atoms_new_colvar_forces = "+cvm::to_str(atoms_new_colvar_forces,
                                                    cvmodule->cv_width,
                                                    cvmodule->cv_prec)+"\n");
  cvmodule->log(cvm::line_marker);

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "atom_groups_new_colvar_forces = "+
           cvm::to_str(atom_groups_new_colvar_forces,
                       cvmodule->cv_width,
                       cvmodule->cv_prec)+"\n");

  cvmodule->log(cvm::line_marker);

  cvmodule->log("Step "+cvm::to_str(cvmodule->step_absolute())+", "+
           "volmaps_new_colvar_forces = "+
           cvm::to_str(volmaps_new_colvar_forces)+"\n");

  cvmodule->log(cvm::line_marker);
}


void colvarproxy::log(std::string const &message)
{
  std::fprintf(stdout, "colvars: %s", message.c_str());
}


void colvarproxy::error(std::string const &message)
{
  // TODO handle errors?
  colvarproxy::log(message);
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


