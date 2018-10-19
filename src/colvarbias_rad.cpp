// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias_rad.h"


colvarbias_rad::colvarbias_rad(char const *key)
  : colvarbias(key)
{
  kernel_coupling_time = 0.0;
  kernel_type = kt_none;

  kernel_num_samples = 0;
}


int colvarbias_rad::init_centers(std::string const &conf)
{
  size_t i;

  bool null_centers = true;
  if (colvar_centers.size() == 0) {
    colvar_centers.resize(colvars.size());
    for (i = 0; i < colvars.size(); i++) {
      colvar_centers[i].type(colvars[i]->value());
      colvar_centers[i].reset();
    }
  }

  if (get_keyval(conf, "centers", colvar_centers, colvar_centers)) {
    for (i = 0; i < colvars.size(); i++) {
      if (cvm::debug()) {
        cvm::log("colvarbias_rad: parsing initial centers, i = "+cvm::to_str(i)+".\n");
      }
      colvar_centers[i].apply_constraints();
    }
    null_centers = false;
  }

  if (null_centers) {
    colvar_centers.clear();
    return cvm::error("Error: must define the centers of the RAD bias.\n",
                      INPUT_ERROR);
  }

  if (colvar_centers.size() != num_variables()) {
    return cvm::error("Error: number of centers does not match "
                      "that of collective variables.\n", INPUT_ERROR);
  }

  return COLVARS_OK;
}


int colvarbias_rad::init(std::string const &conf)
{
  colvarbias::init(conf);

  enable(f_cvb_apply_force);

  if (init_centers(conf) != COLVARS_OK) {
    return cvm::get_error();
  }

  if (colvar_centers_errors.size() == 0) {
    colvar_centers_errors.resize(colvars.size());
    colvar_centers_errors.assign(colvars.size(), 0.0);
  }
  get_keyval(conf, "centersErrors", colvar_centers_errors, colvar_centers_errors);

  std::string kernel_type_str;
  get_keyval(conf, "kernelType", kernel_type_str, to_lower_cppstr(std::string("inverseSqrtTime")));
  kernel_type_str = to_lower_cppstr(kernel_type_str);
  if (kernel_type_str == to_lower_cppstr(std::string("inverseSqrtTime"))) {
    kernel_type = kt_inv_sqrt_time;
  } else if (kernel_type_str == to_lower_cppstr(std::string("uniform"))) {
    kernel_type = kt_uniform;
  } else if (kernel_type_str == to_lower_cppstr(std::string("none"))) {
    kernel_type = kt_none;
  }

  switch (kernel_type) {
  case kt_none:
    cvm::error("Error: undefined kernel type.\n", INPUT_ERROR);
    return INPUT_ERROR;
    break;
  case kt_inv_sqrt_time:
  case kt_uniform:
  case kt_ntot:
    provide(f_cvb_history_dependent);
    break;
  }

  get_keyval(conf, "couplingTime", kernel_coupling_time, 0.0);
  if (kernel_coupling_time > 0.0) {
    enable(f_cvb_history_dependent);
  } else {
    return cvm::error("Error: a positive couplingTime must be provided.\n",
                      INPUT_ERROR);
  }

  // set up the total deviations
  if (colvar_total_deviations.size() == 0) {
    colvar_total_deviations.resize(colvars.size());
    for (size_t i = 0; i < colvars.size(); i++) {
      colvar_total_deviations[i].type(colvars[i]->value());
      colvar_total_deviations[i].is_derivative(); // constraints are not applied
      colvar_total_deviations[i].reset();
    }
  }

  return COLVARS_OK;
}


int colvarbias_rad::update()
{
  bias_energy = 0.0;

  if (cvm::debug())
    cvm::log("Updating the RAD bias \""+this->name+"\".\n");

  cvm::real weight = 1.0;
  switch (kernel_type) {
  case kt_inv_sqrt_time:
    weight = 1.0/std::sqrt(kernel_coupling_time * (kernel_coupling_time + cvm::step_absolute() * cvm::dt()));
    break;
  case kt_uniform:
    weight = 1.0/kernel_coupling_time;
    break;
  case kt_none:
  case kt_ntot:
    break;
  }

  cvm::real const kB_T = cvm::boltzmann() * cvm::temperature();

  size_t i;
  for (i = 0; i < colvars.size(); i++) {

    cvm::real const unit_scale = 1.0/(colvars[i]->width * colvars[i]->width);

    // shift the center to correct for the external error
    colvarvalue const error_drift = unit_scale *
      colvar_total_deviations[i] *
      (colvar_centers_errors[i]*colvar_centers_errors[i]);
    colvarvalue const corrected_center =
      (colvar_centers[i] + error_drift).apply_constraints();
    colvarvalue const deviation =
      0.5 * colvars[i]->dist2_lgrad(colvars[i]->value(), corrected_center);

    colvar_total_deviations[i] += weight * deviation * cvm::dt();
    bias_energy += kB_T * unit_scale * colvar_total_deviations[i] * deviation;
    colvar_forces[i] = -1.0 * kB_T * unit_scale * colvar_total_deviations[i];
  }

  return COLVARS_OK;
}


std::string const colvarbias_rad::get_state_params() const
{
  std::ostringstream os;
  os << "num_samples " << kernel_num_samples << "\n";
  os << "total_deviations ";
  for (size_t i = 0; i < colvars.size(); i++) {
    os << colvar_total_deviations[i];
  }
  os << "\n";
  return (colvarbias::get_state_params() + os.str());
}


int colvarbias_rad::set_state_params(std::string const &state_conf)
{
  int error_code = COLVARS_OK;

  error_code |= colvarbias::set_state_params(state_conf);

  if (!get_keyval(state_conf, "num_samples", kernel_num_samples)) {
    error_code |= cvm::error("Error: num_samples missing from the restart.\n",
                             INPUT_ERROR);
  }

  if (!get_keyval(state_conf, "total_deviations", colvar_total_deviations)) {
    error_code |= cvm::error("Error: total_deviations missing from the restart.\n",
                             INPUT_ERROR);
  }
  return error_code;
}
