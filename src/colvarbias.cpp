// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <fstream>
#include <cstring>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"
#include "colvarbias.h"
#include "colvargrid.h"


colvarbias::colvarbias(char const *key)
{
  bias_type = colvarparse::to_lower_cppstr(key);
  state_keyword = bias_type;

  rank = -1;
  description = "uninitialized " + bias_type + " bias";

  colvarbias::init_dependencies();

  time_step_factor = 1;

  has_data = false;
  b_output_energy = false;
  output_freq = cvm::restart_out_freq;

  colvarbias::reset();
  state_file_step = 0L;
  matching_state = false;
  biasing_force_scaling_factors = NULL;
}


int colvarbias::init(std::string const &conf)
{
  name = bias_type + cvm::to_str(rank);
  colvarparse::set_string(conf);

  size_t i = 0;

  if (num_variables() == 0) {
    // First initialization

    cvm::log("Initializing a new \""+bias_type+"\" instance.\n");

    // Only allow setting a non-default name on first init
    get_keyval(conf, "name", name, name);

    colvarbias *bias_with_name = cvm::bias_by_name(this->name);
    if (bias_with_name != NULL) {
      if ((bias_with_name->rank != this->rank) ||
          (bias_with_name->bias_type != this->bias_type)) {
        cvm::error("Error: this bias cannot have the same name, \""+this->name+
                   "\", as another bias.\n", INPUT_ERROR);
        return INPUT_ERROR;
      }
    }
    description = "bias " + name;

    {
      // lookup the associated colvars
      std::vector<std::string> colvar_names;
      if (get_keyval(conf, "colvars", colvar_names)) {
        if (num_variables()) {
          cvm::error("Error: cannot redefine the colvars that a bias was already defined on.\n",
                     INPUT_ERROR);
          return INPUT_ERROR;
        }
        for (i = 0; i < colvar_names.size(); i++) {
          add_colvar(colvar_names[i]);
        }
      }
    }

    if (!num_variables()) {
      cvm::error("Error: no collective variables specified.\n", INPUT_ERROR);
      return INPUT_ERROR;
    }

  } else {
    cvm::log("Reinitializing bias \""+name+"\".\n");
  }

  colvar_values.resize(num_variables());
  for (i = 0; i < num_variables(); i++) {
    colvar_values[i].type(colvars[i]->value().type());
    colvar_forces[i].type(colvar_values[i].type());
    previous_colvar_forces[i].type(colvar_values[i].type());
  }

  output_prefix = cvm::output_prefix();

  get_keyval_feature(this, conf, "stepZeroData", f_cvb_step_zero_data, is_enabled(f_cvb_step_zero_data));

  // Write energy to traj file?
  get_keyval(conf, "outputEnergy", b_output_energy, b_output_energy);

  // How often to write full output files?
  get_keyval(conf, "outputFreq", output_freq, output_freq);

  // Disabled by default in base class; default value can be overridden by derived class constructor
  get_keyval_feature(this, conf, "bypassExtendedLagrangian", f_cvb_bypass_ext_lagrangian, is_enabled(f_cvb_bypass_ext_lagrangian), parse_echo);

  get_keyval(conf, "timeStepFactor", time_step_factor, time_step_factor);
  if (time_step_factor < 1) {
    cvm::error("Error: timeStepFactor must be 1 or greater.\n");
    return COLVARS_ERROR;
  }

  // Use the scaling factors from a grid?
  get_keyval_feature(this, conf, "scaledBiasingForce",
                     f_cvb_scale_biasing_force,
                     is_enabled(f_cvb_scale_biasing_force), parse_echo);
  if (is_enabled(f_cvb_scale_biasing_force)) {
    std::string biasing_force_scaling_factors_in_filename;
    get_keyval(conf, "scaledBiasingForceFactorsGrid",
               biasing_force_scaling_factors_in_filename, std::string());
    std::ifstream is;
    cvm::log("Reading scaling factors for the forces of bias " +
             name + " from " + biasing_force_scaling_factors_in_filename);
    is.open(biasing_force_scaling_factors_in_filename.c_str());
    if (!is.is_open()) {
      cvm::error("Error opening the grid file " +
      biasing_force_scaling_factors_in_filename + " for reading");
    }
    biasing_force_scaling_factors = new colvar_grid_scalar(colvars);
    biasing_force_scaling_factors->read_multicol(is, true);
    biasing_force_scaling_factors_bin.assign(num_variables(), 0);
    is.close();
  }

  // Now that children are defined, we can solve dependencies
  enable(f_cvb_active);
  if (cvm::debug()) print_state();

  return COLVARS_OK;
}


int colvarbias::init_dependencies() {
  int i;
  if (features().size() == 0) {
    for (i = 0; i < f_cvb_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_cvb_active, "active", f_type_dynamic);
    require_feature_children(f_cvb_active, f_cv_active);

    init_feature(f_cvb_awake, "awake", f_type_static);
    require_feature_self(f_cvb_awake, f_cvb_active);

    init_feature(f_cvb_step_zero_data, "step_zero_data", f_type_user);

    init_feature(f_cvb_apply_force, "apply_force", f_type_user);
    require_feature_children(f_cvb_apply_force, f_cv_gradient);

    init_feature(f_cvb_bypass_ext_lagrangian, "bypass_extended_Lagrangian_coordinates", f_type_user);

    // The exclusion below prevents the inconsistency where biasing forces are applied onto
    // the actual colvar, while total forces are measured on the extended coordinate
    exclude_feature_self(f_cvb_bypass_ext_lagrangian, f_cvb_get_total_force);

    init_feature(f_cvb_get_total_force, "obtain_total_force", f_type_dynamic);
    require_feature_children(f_cvb_get_total_force, f_cv_total_force);

    init_feature(f_cvb_output_acc_work, "output_accumulated_work", f_type_user);
    require_feature_self(f_cvb_output_acc_work, f_cvb_apply_force);

    init_feature(f_cvb_history_dependent, "history_dependent", f_type_static);

    init_feature(f_cvb_time_dependent, "time_dependent", f_type_static);

    init_feature(f_cvb_scalar_variables, "require_scalar_variables", f_type_static);
    require_feature_children(f_cvb_scalar_variables, f_cv_scalar);

    init_feature(f_cvb_calc_pmf, "calculate_a_PMF", f_type_static);

    init_feature(f_cvb_calc_ti_samples, "calculate_TI_samples", f_type_dynamic);
    require_feature_self(f_cvb_calc_ti_samples, f_cvb_get_total_force);
    require_feature_children(f_cvb_calc_ti_samples, f_cv_grid);

    init_feature(f_cvb_write_ti_samples, "write_TI_samples_", f_type_user);
    require_feature_self(f_cvb_write_ti_samples, f_cvb_calc_ti_samples);

    init_feature(f_cvb_write_ti_pmf, "write_TI_PMF", f_type_user);
    require_feature_self(f_cvb_write_ti_pmf, f_cvb_calc_ti_samples);

    init_feature(f_cvb_scale_biasing_force, "scale_biasing_force", f_type_user);
    require_feature_children(f_cvb_scale_biasing_force, f_cv_grid);

    // check that everything is initialized
    for (i = 0; i < colvardeps::f_cvb_ntot; i++) {
      if (is_not_set(i)) {
        cvm::error("Uninitialized feature " + cvm::to_str(i) + " in " + description);
      }
    }
  }

  // Initialize feature_states for each instance
  feature_states.reserve(f_cvb_ntot);
  for (i = 0; i < f_cvb_ntot; i++) {
    feature_states.push_back(feature_state(true, false));
    // Most features are available, so we set them so
    // and list exceptions below
  }

  // only compute TI samples when deriving from colvarbias_ti
  feature_states[f_cvb_calc_ti_samples].available = false;

  // The feature f_cvb_bypass_ext_lagrangian is only implemented by some derived classes
  // (initially, harmonicWalls)
  feature_states[f_cvb_bypass_ext_lagrangian].available = false;
  // disabled by default; can be changed by derived classes that implement it
  feature_states[f_cvb_bypass_ext_lagrangian].enabled = false;

  return COLVARS_OK;
}


int colvarbias::reset()
{
  bias_energy = 0.0;
  for (size_t i = 0; i < num_variables(); i++) {
    colvar_forces[i].reset();
  }
  return COLVARS_OK;
}


colvarbias::colvarbias()
  : colvarparse(), has_data(false)
{}


colvarbias::~colvarbias()
{
  colvarbias::clear();
}


int colvarbias::clear()
{
  free_children_deps();

  // Remove references to this bias from colvars
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       ++cvi) {
    for (std::vector<colvarbias *>::iterator bi = (*cvi)->biases.begin();
         bi != (*cvi)->biases.end();
         ++bi) {
      if ( *bi == this) {
        (*cvi)->biases.erase(bi);
        break;
      }
    }
  }

  colvarmodule *cv = cvm::main();
  // ...and from the colvars module
  for (std::vector<colvarbias *>::iterator bi = cv->biases.begin();
       bi != cv->biases.end();
       ++bi) {
    if ( *bi == this) {
      cv->biases.erase(bi);
      break;
    }
  }

  if (biasing_force_scaling_factors != NULL) {
    delete biasing_force_scaling_factors;
    biasing_force_scaling_factors = NULL;
    biasing_force_scaling_factors_bin.clear();
  }

  cv->config_changed();

  return COLVARS_OK;
}


int colvarbias::clear_state_data()
{
  // no mutable content to delete for base class
  return COLVARS_OK;
}


int colvarbias::add_colvar(std::string const &cv_name)
{
  if (colvar *cv = cvm::colvar_by_name(cv_name)) {

    if (cvm::debug()) {
      cvm::log("Applying this bias to collective variable \""+
               cv->name+"\".\n");
    }

    colvars.push_back(cv);
    cv->biases.push_back(this); // add back-reference to this bias to colvar

    // Add dependency link. All biases need at least the value of each colvar
    // although possibly not at all timesteps
    add_child(cv);

    colvar_forces.push_back(colvarvalue());
    colvar_forces.back().type(cv->value()); // make sure each force is initialized to zero
    colvar_forces.back().is_derivative(); // colvar constraints are not applied to the force
    colvar_forces.back().reset();
    previous_colvar_forces.push_back(colvar_forces.back());

  } else {
    cvm::error("Error: cannot find a colvar named \""+
               cv_name+"\".\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  return COLVARS_OK;
}


int colvarbias::update()
{
  if (cvm::debug()) {
    cvm::log("Updating the "+bias_type+" bias \""+this->name+"\".\n");
  }

  int error_code = COLVARS_OK;

  has_data = true;

  // Update the cached colvar values
  for (size_t i = 0; i < num_variables(); i++) {
    colvar_values[i] = colvars[i]->value();
  }

  error_code |= calc_energy(NULL);
  error_code |= calc_forces(NULL);

  return error_code;
}


bool colvarbias::can_accumulate_data()
{
  colvarproxy *proxy = cvm::main()->proxy;
  if (((cvm::step_relative() > 0) && !proxy->simulation_continuing()) ||
      is_enabled(f_cvb_step_zero_data)) {
    return true;
  }
  return false;
}


int colvarbias::calc_energy(std::vector<colvarvalue> const *)
{
  bias_energy = 0.0;
  return COLVARS_OK;
}


int colvarbias::calc_forces(std::vector<colvarvalue> const *)
{
  for (size_t ir = 0; ir < num_variables(); ir++) {
    colvar_forces[ir].reset();
  }
  return COLVARS_OK;
}


void colvarbias::communicate_forces()
{
  if (! is_enabled(f_cvb_apply_force)) {
    return;
  }
  cvm::real biasing_force_factor = 1.0;
  size_t i = 0;
  if (is_enabled(f_cvb_scale_biasing_force)) {
    for (i = 0; i < num_variables(); i++) {
      biasing_force_scaling_factors_bin[i] = biasing_force_scaling_factors->current_bin_scalar(i);
    }
    if (biasing_force_scaling_factors->index_ok(biasing_force_scaling_factors_bin)) {
      biasing_force_factor *= biasing_force_scaling_factors->value(biasing_force_scaling_factors_bin);
    }
  }
  for (i = 0; i < num_variables(); i++) {
    if (cvm::debug()) {
      cvm::log("Communicating a force to colvar \""+
               variables(i)->name+"\".\n");
    }
    // Impulse-style multiple timestep
    // Note that biases with different values of time_step_factor
    // may send forces to the same colvar
    // which is why rescaling has to happen now: the colvar is not
    // aware of this bias' time_step_factor
    if (is_enabled(f_cvb_bypass_ext_lagrangian)) {
      variables(i)->add_bias_force_actual_value(cvm::real(time_step_factor) * colvar_forces[i] * biasing_force_factor);
    } else {
      variables(i)->add_bias_force(cvm::real(time_step_factor) * colvar_forces[i] * biasing_force_factor);
    }
    previous_colvar_forces[i] = colvar_forces[i];
  }
}


int colvarbias::end_of_step()
{
  return COLVARS_OK;
}


int colvarbias::change_configuration(std::string const & /* conf */)
{
  cvm::error("Error: change_configuration() not implemented.\n",
             COLVARS_NOT_IMPLEMENTED);
  return COLVARS_NOT_IMPLEMENTED;
}


cvm::real colvarbias::energy_difference(std::string const & /* conf */)
{
  cvm::error("Error: energy_difference() not implemented.\n",
             COLVARS_NOT_IMPLEMENTED);
  return 0.0;
}


// So far, these are only implemented in colvarbias_abf
int colvarbias::bin_num()
{
  cvm::error("Error: bin_num() not implemented.\n");
  return COLVARS_NOT_IMPLEMENTED;
}
int colvarbias::current_bin()
{
  cvm::error("Error: current_bin() not implemented.\n");
  return COLVARS_NOT_IMPLEMENTED;
}
int colvarbias::bin_count(int /* bin_index */)
{
  cvm::error("Error: bin_count() not implemented.\n");
  return COLVARS_NOT_IMPLEMENTED;
}
int colvarbias::replica_share()
{
  cvm::error("Error: replica_share() not implemented.\n");
  return COLVARS_NOT_IMPLEMENTED;
}


std::string const colvarbias::get_state_params() const
{
  std::ostringstream os;
  os << "step " << cvm::step_absolute() << "\n"
     << "name " << this->name << "\n";
  return os.str();
}


int colvarbias::set_state_params(std::string const &conf)
{
  matching_state = false;

  std::string check_name = "";
  colvarparse::get_keyval(conf, "name", check_name,
                          std::string(""), colvarparse::parse_silent);

  if (check_name.size() == 0) {
    cvm::error("Error: \""+bias_type+"\" block within the restart file "
               "has no identifiers.\n", INPUT_ERROR);
  }

  if (check_name != this->name) {
    if (cvm::debug()) {
      cvm::log("Ignoring state of bias \""+check_name+
               "\": this bias is named \""+name+"\".\n");
    }
    return COLVARS_OK;
  }

  matching_state = true;

  colvarparse::get_keyval(conf, "step", state_file_step,
                          cvm::step_absolute(), colvarparse::parse_silent);

  return COLVARS_OK;
}


std::ostream & colvarbias::write_state(std::ostream &os)
{
  if (cvm::debug()) {
    cvm::log("Writing state file for bias \""+name+"\"\n");
  }
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(cvm::cv_prec);
  os << state_keyword << " {\n"
     << "  configuration {\n";
  std::istringstream is(get_state_params());
  std::string line;
  while (std::getline(is, line)) {
    os << "    " << line << "\n";
  }
  os << "  }\n";
  write_state_data(os);
  os << "}\n\n";
  return os;
}


std::istream & colvarbias::read_state(std::istream &is)
{
  std::streampos const start_pos = is.tellg();

  std::string key, brace, conf;
  if ( !(is >> key)   || !(key == state_keyword || key == bias_type) ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block("configuration", &conf)) ||
       (set_state_params(conf) != COLVARS_OK) ) {
    cvm::error("Error: in reading state configuration for \""+bias_type+
               "\" bias \""+
               this->name+"\" at position "+
               cvm::to_str(static_cast<size_t>(is.tellg()))+
               " in stream.\n", INPUT_ERROR);
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  if (matching_state == false) {
    // This state is not for this bias
    is.seekg(start_pos, std::ios::beg);
    return is;
  }

  cvm::log("Restarting "+bias_type+" bias \""+name+"\" from step number "+
           cvm::to_str(state_file_step)+".\n");

  if (!read_state_data(is)) {
    cvm::error("Error: in reading state data for \""+bias_type+"\" bias \""+
               this->name+"\" at position "+
               cvm::to_str(static_cast<size_t>(is.tellg()))+
               " in stream.\n", INPUT_ERROR);
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
  }

  is >> brace;
  if (brace != "}") {
    cvm::error("Error: corrupt restart information for \""+bias_type+"\" bias \""+
               this->name+"\": no matching brace at position "+
               cvm::to_str(static_cast<size_t>(is.tellg()))+
               " in stream.\n");
    is.setstate(std::ios::failbit);
  }

  return is;
}


int colvarbias::write_state_prefix(std::string const &prefix)
{
  std::string const filename =
    cvm::state_file_prefix(prefix.c_str())+".colvars.state";
  std::ostream *os = cvm::proxy->output_stream(filename.c_str());
  int error_code = COLVARS_OK;
  if (os != NULL) {
    os->setf(std::ios::scientific, std::ios::floatfield);
    error_code = write_state(*os).good() ? COLVARS_OK : FILE_ERROR;
  } else {
    error_code = FILE_ERROR;
  }
  cvm::proxy->close_output_stream(filename.c_str());
  return error_code;
}


int colvarbias::write_state_string(std::string &output)
{
  std::ostringstream os;
  if (!write_state(os)) {
    return cvm::error("Error: in writing state of bias \""+name+
                      "\" to buffer.\n", FILE_ERROR);
  }
  output = os.str();
  return COLVARS_OK;
}


int colvarbias::read_state_prefix(std::string const &prefix)
{
  std::string filename((prefix+std::string(".colvars.state")).c_str());
  std::ifstream is(filename.c_str());
  if (!is.good()) {
    // try without the suffix
    is.clear();
    filename = prefix;
    is.open(filename.c_str());
  }
  return read_state(is).good() ? COLVARS_OK :
    cvm::error("Error: in reading state for \""+name+"\" from input file \""+
               std::string(filename)+"\".\n", FILE_ERROR);
}


int colvarbias::read_state_string(char const *buffer)
{
  if (buffer != NULL) {
    size_t const buffer_size = strlen(buffer);
    if (cvm::debug()) {
      cvm::log("colvarbias::read_state_string() with argument:\n");
      cvm::log(buffer);
    }

    if (buffer_size > 0) {
      std::istringstream is;
      is.rdbuf()->pubsetbuf(const_cast<char *>(buffer), buffer_size);
      return read_state(is).good() ? COLVARS_OK :
        cvm::error("Error: in reading state for \""+name+"\" from buffer.\n",
                   FILE_ERROR);
    }
    return COLVARS_OK;
  }
  return cvm::error("Error: NULL pointer for colvarbias::read_state_string()",
                    BUG_ERROR);
}


std::istream & colvarbias::read_state_data_key(std::istream &is, char const *key)
{
  std::streampos const start_pos = is.tellg();
  std::string key_in;
  if ( !(is >> key_in) ||
       !(to_lower_cppstr(key_in) == to_lower_cppstr(std::string(key))) ) {
    cvm::error("Error: in reading restart configuration for "+
               bias_type+" bias \""+this->name+"\" at position "+
               cvm::to_str(static_cast<size_t>(is.tellg()))+
               " in stream.\n", INPUT_ERROR);
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }
  return is;
}



std::ostream & colvarbias::write_traj_label(std::ostream &os)
{
  os << " ";
  if (b_output_energy)
    os << " E_"
       << cvm::wrap_string(this->name, cvm::en_width-2);
  return os;
}


std::ostream & colvarbias::write_traj(std::ostream &os)
{
  os << " ";
  if (b_output_energy)
    os << " "
       << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width)
       << bias_energy;
  return os;
}



colvarbias_ti::colvarbias_ti(char const *key)
  : colvarbias(key)
{
  colvarproxy *proxy = cvm::main()->proxy;
  provide(f_cvb_calc_ti_samples);
  if (!proxy->total_forces_same_step()) {
    // Samples at step zero can not be collected
    feature_states[f_cvb_step_zero_data].available = false;
  }
  ti_avg_forces = NULL;
  ti_count = NULL;
}


colvarbias_ti::~colvarbias_ti()
{
  colvarbias_ti::clear_state_data();
}


int colvarbias_ti::clear_state_data()
{
  if (ti_avg_forces != NULL) {
    delete ti_avg_forces;
    ti_avg_forces = NULL;
  }
  if (ti_count != NULL) {
    delete ti_count;
    ti_count = NULL;
  }
  return COLVARS_OK;
}


int colvarbias_ti::init(std::string const &conf)
{
  int error_code = COLVARS_OK;

  get_keyval_feature(this, conf, "writeTISamples",
                     f_cvb_write_ti_samples,
                     is_enabled(f_cvb_write_ti_samples));

  get_keyval_feature(this, conf, "writeTIPMF",
                     f_cvb_write_ti_pmf,
                     is_enabled(f_cvb_write_ti_pmf));

  if ((num_variables() > 1) && is_enabled(f_cvb_write_ti_pmf)) {
    return cvm::error("Error: only 1-dimensional PMFs can be written "
                      "on the fly.\n"
                      "Consider using writeTISamples instead and "
                      "post-processing the sampled free-energy gradients.\n",
                      COLVARS_NOT_IMPLEMENTED);
  } else {
    error_code |= init_grids();
  }

  if (is_enabled(f_cvb_write_ti_pmf)) {
    enable(f_cvb_write_ti_samples);
  }

  if (is_enabled(f_cvb_calc_ti_samples)) {
    std::vector<std::string> const time_biases =
      cvm::main()->time_dependent_biases();
    if (time_biases.size() > 0) {
      if ((time_biases.size() > 1) || (time_biases[0] != this->name)) {
        for (size_t i = 0; i < num_variables(); i++) {
          if (! variables(i)->is_enabled(f_cv_subtract_applied_force)) {
            return cvm::error("Error: cannot collect TI samples while other "
                              "time-dependent biases are active and not all "
                              "variables have subtractAppliedForces on.\n",
                              INPUT_ERROR);
          }
        }
      }
    }
  }

  if (is_enabled(f_cvb_write_ti_pmf) || is_enabled(f_cvb_write_ti_samples)) {
    cvm::main()->cite_feature("Internal-forces free energy estimator");
  }

  return error_code;
}


int colvarbias_ti::init_grids()
{
  if (is_enabled(f_cvb_calc_ti_samples)) {
    if (ti_avg_forces == NULL) {
      ti_bin.resize(num_variables());
      ti_system_forces.resize(num_variables());
      for (size_t icv = 0; icv < num_variables(); icv++) {
        ti_system_forces[icv].type(variables(icv)->value());
        ti_system_forces[icv].is_derivative();
        ti_system_forces[icv].reset();
      }
      ti_avg_forces = new colvar_grid_gradient(colvars);
      ti_count = new colvar_grid_count(colvars);
      ti_avg_forces->samples = ti_count;
      ti_count->has_parent_data = true;
    }
  }

  return COLVARS_OK;
}


int colvarbias_ti::update()
{
  return update_system_forces(NULL);
}


int colvarbias_ti::update_system_forces(std::vector<colvarvalue> const
                                        *subtract_forces)
{
  if (! is_enabled(f_cvb_calc_ti_samples)) {
    return COLVARS_OK;
  }

  has_data = true;

  if (cvm::debug()) {
    cvm::log("Updating system forces for bias "+this->name+"\n");
  }

  colvarproxy *proxy = cvm::main()->proxy;

  size_t i;

  if (proxy->total_forces_same_step()) {
    for (i = 0; i < num_variables(); i++) {
      ti_bin[i] = ti_avg_forces->current_bin_scalar(i);
    }
  }

  // Collect total colvar forces
  if ((cvm::step_relative() > 0) || proxy->total_forces_same_step()) {
    if (ti_avg_forces->index_ok(ti_bin)) {
      for (i = 0; i < num_variables(); i++) {
        if (variables(i)->is_enabled(f_cv_subtract_applied_force)) {
          // this colvar is already subtracting all applied forces
          ti_system_forces[i] = variables(i)->total_force();
        } else {
          ti_system_forces[i] = variables(i)->total_force() -
            ((subtract_forces != NULL) ?
             (*subtract_forces)[i] : previous_colvar_forces[i]);
        }
      }
      if (cvm::step_relative() > 0 || is_enabled(f_cvb_step_zero_data)) {
        ti_avg_forces->acc_value(ti_bin, ti_system_forces);
      }
    }
  }

  if (!proxy->total_forces_same_step()) {
    // Set the index for use in the next iteration, when total forces come in
    for (i = 0; i < num_variables(); i++) {
      ti_bin[i] = ti_avg_forces->current_bin_scalar(i);
    }
  }

  return COLVARS_OK;
}


std::string const colvarbias_ti::get_state_params() const
{
  return std::string("");
}


int colvarbias_ti::set_state_params(std::string const & /* state_conf */)
{
  return COLVARS_OK;
}


std::ostream & colvarbias_ti::write_state_data(std::ostream &os)
{
  if (! is_enabled(f_cvb_calc_ti_samples)) {
    return os;
  }
  os << "\nhistogram\n";
  ti_count->write_raw(os);
  os << "\nsystem_forces\n";
  ti_avg_forces->write_raw(os);
  return os;
}


std::istream & colvarbias_ti::read_state_data(std::istream &is)
{
  if (! is_enabled(f_cvb_calc_ti_samples)) {
    return is;
  }
  if (cvm::debug()) {
    cvm::log("Reading state data for the TI estimator.\n");
  }
  if (! read_state_data_key(is, "histogram")) {
    return is;
  }
  if (! ti_count->read_raw(is)) {
    return is;
  }
  if (! read_state_data_key(is, "system_forces")) {
    return is;
  }
  if (! ti_avg_forces->read_raw(is)) {
    return is;
  }
  if (cvm::debug()) {
    cvm::log("Done reading state data for the TI estimator.\n");
  }
  return is;
}


int colvarbias_ti::write_output_files()
{
  if (!has_data) {
    // nothing to write
    return COLVARS_OK;
  }

  std::string const ti_output_prefix = cvm::output_prefix()+"."+this->name;

  std::ostream *os = NULL;

  if (is_enabled(f_cvb_write_ti_samples)) {
    std::string const ti_count_file_name(ti_output_prefix+".ti.count");
    os = cvm::proxy->output_stream(ti_count_file_name);
    if (os) {
      ti_count->write_multicol(*os);
      cvm::proxy->close_output_stream(ti_count_file_name);
    }

    std::string const ti_grad_file_name(ti_output_prefix+".ti.force");
    os = cvm::proxy->output_stream(ti_grad_file_name);
    if (os) {
      ti_avg_forces->write_multicol(*os);
      cvm::proxy->close_output_stream(ti_grad_file_name);
    }
  }

  if (is_enabled(f_cvb_write_ti_pmf)) {
    std::string const pmf_file_name(ti_output_prefix+".ti.pmf");
    cvm::log("Writing TI PMF to file \""+pmf_file_name+"\".\n");
    os = cvm::proxy->output_stream(pmf_file_name);
    if (os) {
      // get the FE gradient
      ti_avg_forces->multiply_constant(-1.0);
      ti_avg_forces->write_1D_integral(*os);
      ti_avg_forces->multiply_constant(-1.0);
      cvm::proxy->close_output_stream(pmf_file_name);
    }
  }

  return COLVARS_OK;
}


// Static members

std::vector<colvardeps::feature *> colvarbias::cvb_features;
