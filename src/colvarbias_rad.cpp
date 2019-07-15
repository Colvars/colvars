// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"
#include "colvarbias_rad.h"


colvarbias_rad::colvarbias_rad(char const *key)
  : colvarbias(key)
{
  kernel_coupling_time = 0.0;
  kernel_type = kt_none;
  opt_type = opt_none;
  use_norm_1 = false;

  colvar_aver_deviation = 0.0;
  colvar_cum_error = 0.0;

  colvar_rad_steps = 0; // XXX to be deleted

  provide(f_cvb_opt_cv_params);
}


int colvarbias_rad::init(std::string const &conf)
{
  colvarbias::init(conf);

  enable(f_cvb_apply_force);
  size_t i;

  if (init_centers(conf) != COLVARS_OK) {
    return cvm::get_error();
  }

  // colvar_widths_c

  if (colvar_widths_c.size() != num_variables()) {
    colvar_widths_c.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      colvar_widths_c[i] = variables(i)->width;
    }
  }

  get_keyval(conf, "colvarWidths", colvar_widths_c, colvar_widths_c);

  for (i = 0; i < num_variables(); i++) {
     // revert width parameters scaling if required
     colvar_widths_c[i] = variables(i)->rescale_width(colvar_widths_c[i]);
     colvar_widths_c[i] = colvar_widths_c[i]*colvar_widths_c[i];
  }

  get_keyval_feature(this, conf, "optParams",
                     f_cvb_opt_cv_params, is_enabled(f_cvb_opt_cv_params));

  if (is_enabled(f_cvb_opt_cv_params)) {
    cvm::log("NOTE: parameters optimization assumes linear parameters functions or linear approximations \n");
    std::string opt_type_str;
    get_keyval(conf, "ParOptType", opt_type_str, std::string("chisquare"));
    opt_type_str = to_lower_cppstr(opt_type_str);
    if (opt_type_str == to_lower_cppstr(std::string("chisquare"))) {
      opt_type = opt_chisquare;
    } else if (opt_type_str == to_lower_cppstr(std::string("lambda"))) {
      opt_type = opt_lambda;
    } else if (opt_type_str == to_lower_cppstr(std::string("none"))) {
      opt_type = opt_none;
    }

    switch (opt_type) {
    case opt_none:
      cvm::error("Error: undefined parameters optimization type.\n", INPUT_ERROR);
      return INPUT_ERROR;
      break;
    case opt_chisquare:
      // set up the total chideviations
      get_keyval(conf, "optParamsSteps",colvar_chisquare_opt_steps, 1000);
      get_keyval(conf, "optParamsToll",colvar_chisquare_opt_toll, 0.00000000001);
      if (colvar_total_chideviations.size() == 0) {
        colvar_total_chideviations.resize(num_variables());
        for (i = 0; i < num_variables(); i++) {
          colvar_total_chideviations[i].type(variables(i));
          colvar_total_chideviations[i].is_derivative(); // constraints are not applied
          colvar_total_chideviations[i].r000et();
        }
      }
      break;
    case opt_lambda:
      get_keyval(conf, "paramsCouplingTime", params_coupling_time, 1.0);
      break;
    }
  }

  get_keyval(conf, "radOutFreq", rad_out_freq, 1000);

  size_t ii;
  size_t t;
  cvm::real eback;

//  if (is_enabled(f_cvb_opt_cv_params)) {
//    // initialize and get number of parameters and values of parameters for each set of CVs

//    if (val_params.size()==0) {
//      val_params.resize(colvars.size());
//    }

//    // read value of the parameters

//    for (i = 0; i < num_variables(); i++) {
//       variables(i)->get_params(val_params[i]);
//    }
//  }

  if (colvar_deviation.size() == 0) {
    colvar_deviation.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      colvar_deviation[i].type(variables(i)->value());
      colvar_deviation[i].is_derivative(); // constraints are not applied
      colvar_deviation[i].reset();
    }
  }

  get_keyval(conf, "setDevFromExpToOne", fix_chi_square_one, false);

  colvar_errors_scale=1.; // scale factor of the experimental error

  get_keyval(conf, "useNorm1", use_norm_1, use_norm_1);

  if (use_norm_1) {
    cvm::log("Using Norm of type 1 to calculate the deviation from the experimental value\n");
  }

  std::string kernel_type_str;
  get_keyval(conf, "kernelType", kernel_type_str, std::string("inverseSqrtTime"));
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
    colvar_total_deviations.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      colvar_total_deviations[i].type(variables(i));
      colvar_total_deviations[i].is_derivative(); // constraints are not applied
      colvar_total_deviations[i].reset();
    }
  }

  return COLVARS_OK;
}


int colvarbias_rad::init_centers(std::string const &conf)
{
  size_t i;

  if (colvar_centers.size() == 0) {
    colvar_centers.resize(num_variables());
    colvar_exp_centers.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      colvar_centers[i].type(variables(i)->value());
      colvar_centers[i].reset();
      colvar_exp_centers[i].type(variables(i)->value());
      colvar_exp_centers[i].reset();
    }
  }

  // get centers from variables

  for (i = 0; i < num_variables(); i++) {
     variables(i)->get_exp_val(colvar_exp_centers[i]);
     colvar_centers[i]=colvar_exp_centers[i];
  }

  // read centers from input

  get_keyval(conf, "expCenters", colvar_exp_centers, colvar_exp_centers);

  for (i = 0; i < num_variables(); i++) {
     colvar_centers[i]=colvar_exp_centers[i]; // assign initial centers (or optimal values)
  }

  for (i = 0; i < num_variables(); i++) {
    if (cvm::debug()) {
      cvm::log("colvarbias_restraint: parsing initial centers, i = "+cvm::to_str(i)+".\n");
    }
    colvar_centers[i].apply_constraints();
  }

  if (colvar_centers.size() != num_variables()) {
    cvm::error("Error: number of centers does not match "
               "that of collective variables.\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  if (colvar_centers_errors.size() == 0) {
    colvar_centers_errors.resize(colvars.size());
    colvar_centers_errors.assign(colvars.size(), 0.0);
  }
  get_keyval(conf, "centersErrors", colvar_centers_errors, colvar_centers_errors);

  // calc average error

  colvar_cum_error=0;
  for (i = 0; i < num_variables(); i++) {
     colvar_cum_error=colvar_cum_error+colvar_centers_errors[i]*colvar_centers_errors[i];
  }
  colvar_cum_error=colvar_cum_error/double(num_variables());

  return COLVARS_OK;
}


int colvarbias_rad::clear()
{
  int error_code = COLVARS_OK;
  if (rad_out_os) {
    cvm::proxy->close_output_stream(rad_out_file_name());
  }
  if (rad_param_os) {
    cvm::proxy->close_output_stream(rad_param_file_name());
  }
  return error_code | colvarbias::clear();
}



int colvarbias_rad::update()
{
  colvarbias::update();

  if (cvm::debug())
    cvm::log("Updating the linear optimizer bias \""+this->name+"\".\n");

  cvm::real const kBT = cvm::boltzmann() * cvm::temperature();

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

  size_t i;
  size_t t;
  size_t ii;

  colvar_aver_deviation = 0.0;
  cvm::real lambdasum=0.0;
  cvm::real lambda2;
  cvm::real lambda2sum=0.0;
  cvm::real error_fact;

  colvar_rad_steps = colvar_rad_steps+1; // XXX to be deleted

  for (i = 0; i < num_variables(); i++) {
    cvm::real const unit_scale = 1.0/colvar_widths_c[i];
    error_fact = colvar_centers_errors[i]*colvar_centers_errors[i]/colvar_errors_scale;
    colvarvalue const deviation =
      0.5 * variables(i)->dist2_lgrad(variables(i)->value(), colvar_centers[i]);

    colvar_total_deviations[i] += weight * variables(i)->paramscale(deviation) * cvm::dt();
    // scaling removed to calculate the actual force on the CV 
    colvar_forces[i] = -1.0 * kBT * variables(i)->paramscale(unit_scale * colvar_total_deviations[i]);
    bias_energy += -colvar_forces[i]*deviation;

    colvarvalue const error_drift = -(colvar_forces[i]/kBT)*error_fact;
    colvar_centers[i] = colvar_exp_centers[i] + error_drift;
    colvar_centers[i].apply_constraints();
    colvar_deviation[i]=error_drift/colvar_centers_errors[i];
    for (t=0;t<colvar_centers[i].size();t++) {
       lambda2=colvar_forces[i].vector1d_value[t]/kBT;
       lambda2=lambda2*lambda2;
       lambdasum+=colvar_centers_errors[i]*std::sqrt(lambda2);
       lambda2sum+=colvar_centers_errors[i]*colvar_centers_errors[i]*lambda2;;
    }

    switch (opt_type) {
    case opt_chisquare:
      // get running average for optimizing the parameters by minimizing the chisquare
      colvarvalue const chideviation =
        0.5 * variables(i)->dist2_lgrad(variables(i)->value(), colvar_exp_centers[i]);
      colvar_total_chideviations[i] += variables(i)->paramscale(chideviation);
      variables(i)->update_params_rad_chis(colvar_total_chideviations[i]/colvar_rad_steps, colvar_exp_centers[i],
                                       colvar_chisquare_opt_steps,colvar_chisquare_opt_toll);
    case opt_lambda:
      variables(i)->update_params_rad(colvar_forces[i]/kBT, colvar_centers[i],
                                       params_coupling_time, weight, unit_scale,
                                       colvar_widths_c[i]);
    case opt_none:
      break;
    }
  }
  
  if (use_norm_1) {
    for (i = 0; i < num_variables(); i++) {
      cvm::real colvar_local_deviation=0.;
      colvar_deviation[i].set_absolute_value();
      colvarvalue const save_colvar_value = colvar_deviation[i];
      colvar_deviation[i].set_ones();
      colvar_local_deviation=colvar_deviation[i]*save_colvar_value;
      colvar_aver_deviation+=colvar_local_deviation;
    }
  } else {
    for (i = 0; i < num_variables(); i++) {
      cvm::real colvar_local_deviation=0.;
      colvar_local_deviation=colvar_deviation[i]*colvar_deviation[i];
      colvar_aver_deviation+=colvar_local_deviation;
    }
  }


  colvar_aver_deviation=colvar_aver_deviation/double(variables_num_dimensions());

  lambdasum=lambdasum/double(variables_num_dimensions());
  lambda2sum=std::sqrt(lambda2sum/double(variables_num_dimensions()));

  if (!use_norm_1) {
    colvar_aver_deviation=std::sqrt(colvar_aver_deviation);
  }

  if (fix_chi_square_one) {

    colvar_errors_scale=lambda2sum;
    if (use_norm_1) colvar_errors_scale=lambdasum;

  }
  // check

  return COLVARS_OK;
}


int colvarbias_rad::setup_output()
{
  int error_code = COLVARS_OK;

  error_code |= colvarbias::setup_output();

  if (rad_out_os != NULL) {
    return error_code;
  }

  if (rad_out_freq) {

    cvm::proxy->backup_file(rad_out_file_name());
    rad_out_os = cvm::proxy->output_stream(rad_out_file_name());
    if (!rad_out_os) return FILE_ERROR;

    if (is_enabled(f_cvb_opt_cv_params)) {
      cvm::proxy->backup_file(rad_param_file_name());
      rad_param_os = cvm::proxy->output_stream(rad_param_file_name());
      if (!rad_param_os) return FILE_ERROR;
    }
  }
}


int colvarbias_rad::write_traj_files()
{
  int error_code = COLVARS_OK;

  error_code |= colvarbias::write_traj_files();

  if ((cvm::step_absolute() % rad_out_freq) == 0) {
    if (rad_out_os) {
      std::ostream &os = *rad_out_os;
      os.setf(std::ios::scientific, std::ios::floatfield);
      os << std::setw(cvm::it_width)
         << cvm::step_absolute()
         << "  "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_aver_deviation
         << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << std::sqrt(colvar_cum_error/colvar_errors_scale)
         << " "
         << std::setprecision(cvm::it_width)
         << variables_num_dimensions()
         << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_cum_error
         << "\n";
      error_code |= (os.good() ? COLVARS_OK : FILE_ERROR);
    }

    if (rad_param_os) {
      std::ostream &os = *rad_param_os;
      os.setf(std::ios::scientific, std::ios::floatfield);
      for (int i = 0; i < num_variables(); i++) {
         os << std::setw(cvm::it_width) << cvm::step_absolute()
            << cvm::step_absolute()
            << " ";
         os << variables(i)->name
            << " ";
         variables(i)->write_params_rad(os);
      }
      os << "\n";
      error_code |= (os.good() ? COLVARS_OK : FILE_ERROR);
    }
  }

  return error_code;
}


std::string const colvarbias_rad::get_state_params() const
{
  std::ostringstream os;
  os << "total_deviations ";
  for (size_t i = 0; i < num_variables(); i++) {
    os << colvar_total_deviations[i];
  }
  if (fix_chi_square_one) {
    os << "    error_scale " << colvar_errors_scale << "\n";
  }

  os << "\n";
  return (colvarbias::get_state_params() + os.str());
}


int colvarbias_rad::set_state_params(std::string const &state_conf)
{
  int error_code = COLVARS_OK;

  error_code |= colvarbias::set_state_params(state_conf);

  if (!get_keyval(state_conf, "total_deviations", colvar_total_deviations)) {
    error_code |= cvm::error("Error: total_deviations missing from the restart.\n",
                             INPUT_ERROR);
  }

  if (fix_chi_square_one) {
    if (!get_keyval(state_conf, "error_scale", colvar_errors_scale))
      cvm::error("Error: error_scale missing from the restart.\n");
  }

  return error_code;
}
