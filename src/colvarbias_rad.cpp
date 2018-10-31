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

  colvar_aver_deviation = 0.0;
  colvar_cum_error = 0.0;
  colvar_cum_uscale = 0.0;
  colvar_maxerror_scale = 0.0;
}


int colvarbias_rad::init(std::string const &conf)
{
  colvarbias::init(conf);

  enable(f_cvb_apply_force);
  size_t i;

  if (init_centers(conf) != COLVARS_OK) {
    return cvm::get_error();
  }

  cvm::real colvar_single_e;
  if (get_keyval(conf, "oneCenterError", colvar_single_e, colvar_single_e)) {
    colvar_centers_errors.assign(colvars.size(), colvar_single_e);
  }

  // calc average error

  colvar_cum_error=0;
  for (i = 0; i < colvars.size(); i++) {
     colvar_cum_error=colvar_cum_error+colvar_centers_errors[i]*colvar_centers_errors[i];
  }
  colvar_cum_error=colvar_cum_error/double(colvars.size());

  // colvar_widths_c

  if (colvar_widths_c.size() == 0) {
    colvar_widths_c.resize(colvars.size());
    colvar_widths_c.assign(colvars.size(), 1.0);
  }

  get_keyval(conf, "colvarWidths", colvar_widths_c, colvar_widths_c);

  cvm::real colvar_single_w;

  if (get_keyval(conf, "oneColvarWidth", colvar_single_w, colvar_single_w)) {
    colvar_widths_c.assign(colvars.size(), colvar_single_w);
  }

  colvar_cum_uscale=0;
  colvar_size_tot=0;
  for (i = 0; i < colvars.size(); i++) {
     colvar_cum_uscale=colvar_cum_uscale+colvar_widths_c[i]*colvar_widths_c[i];
     colvar_size_tot+=colvar_centers[i].size();
  }

  colvar_cum_uscale=colvar_cum_uscale/double(colvars.size());
  colvar_cum_uscale=1/colvar_cum_uscale;

  // set colvar types

  if (colvar_types.size() == 0) {
    colvar_types.resize(colvars.size());
    for (i = 0; i < colvars.size(); i++) {
       colvar_types[i]=to_lower_cppstr(std::string("generic"));
    }
  }

  int num_types_cvs=2; // just deer variables and generic for now

  if (numtypes.size() == 0) {
    numtypes.resize(num_types_cvs);
    for (i = 0; i < num_types_cvs; i++) {
       numtypes[i]=0;
    }
  }

  whichtypes.resize(num_types_cvs,colvars.size());

  get_keyval(conf, "colvarTypes", colvar_types, colvar_types);
  for (i = 0; i < colvars.size(); i++) {
     if (colvar_types[i]==to_lower_cppstr(std::string("generic"))) {
       numtypes[0]++;
       whichtypes[0][numtypes[0]]=i;
     }
     if (colvar_types[i]==to_lower_cppstr(std::string("deer"))) {
       numtypes[1]++;
       whichtypes[1][numtypes[1]]=i;
     }
  }

  // set colvar types done

  get_keyval(conf, "optParams", opt_params, false);

  get_keyval(conf, "radOutFreq", rad_out_freq, 1000);

  size_t ii;
  size_t t;
  cvm::real eback;
  if (opt_params) {
    // read initial value of the parameters
    if (numtypes[1]>0) {
      mdepth_deer.resize(numtypes[1]);
      alpha_deer.resize(numtypes[1]);
      for (i = 0; i < numtypes[1]; i++) {
         mdepth_deer[i]=0.02; //default value
         alpha_deer[i]=0.0001; //default value
      }
      get_keyval(conf, "deerMdepth", mdepth_deer, mdepth_deer);
      get_keyval(conf, "deerBackAlpha", alpha_deer, alpha_deer);
      for (i = 0; i < numtypes[1]; i++) {
         ii=whichtypes[1][i];
         for (t=0;t<colvar_centers[ii].size();t++) {
            eback=exp(-alpha_deer[i]*sqrt(colvar_times[ii].vector1d_value[t]*colvar_times[ii].vector1d_value[t]));
            colvar_centers[ii].vector1d_value[t]=(colvar_expval[ii].vector1d_value[t]-(1.0-mdepth_deer[i])*eback)/(eback*mdepth_deer[i]);
         }
         colvar_orig_centers[ii]=colvar_centers[ii];
      }
    }
    // to add other variables with parameters to be optimized
    // if (numtypes[2]>0) {
    // ....
    // }
  }

  if (colvar_deviation.size() == 0) {
    colvar_deviation.resize(colvars.size());
    for (i = 0; i < colvars.size(); i++) {
      colvar_deviation[i].type(colvars[i]->value());
      colvar_deviation[i].is_derivative(); // constraints are not applied
      colvar_deviation[i].reset();
    }
  }

  get_keyval(conf, "fixChiSquareOne", fix_chi_square_one, false);

  get_keyval(conf, "maxErrorScale", colvar_maxerror_scale, 100000.);

  colvar_errors_scale=1.; // scale factor of the experimental error

  get_keyval(conf, "useNorm1", use_norm_1, false);

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
    cvm::error("Error: a positive couplingTime must be provided.\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  get_keyval(conf, "paramsCouplingTime", params_coupling_time, 100.0);

  // set up the total deviations
  if (colvar_total_deviations.size() == 0) {
    colvar_total_deviations.resize(colvars.size());
    for (i = 0; i < colvars.size(); i++) {
      colvar_total_deviations[i].type(colvars[i]->value());
      colvar_total_deviations[i].is_derivative(); // constraints are not applied
      colvar_total_deviations[i].reset();
    }
  }

  return COLVARS_OK;
}


int colvarbias_rad::init_centers(std::string const &conf)
{
  size_t i;

  bool null_centers = true;
  bool read_exp_file = false;

  if (colvar_centers.size() == 0) {
    colvar_centers.resize(colvars.size());
    colvar_orig_centers.resize(colvars.size());
    for (i = 0; i < colvars.size(); i++) {
      colvar_centers[i].type(colvars[i]->value());
      colvar_centers[i].reset();
      colvar_orig_centers[i].type(colvars[i]->value());
      colvar_orig_centers[i].reset();
    }
  }

  if (time_files.size() == 0) {
    time_files.resize(colvars.size());
  }

  if (get_keyval(conf, "expFiles", time_files)) {
    read_exp_file = true;
    cvm::log("Initial centers read from expFiles).\n");
    if (colvar_times.size() == 0) {
      colvar_times.resize(colvars.size());
      for (i = 0; i < colvars.size(); i++) {
        colvar_times[i].type(colvars[i]->value());
        colvar_times[i].is_derivative(); // constraints are not applied
        colvar_times[i].reset();
      }
    }

    if (colvar_expval.size() == 0) {
      colvar_expval.resize(colvars.size());
      for (i = 0; i < colvars.size(); i++) {
        colvar_expval[i].type(colvars[i]->value());
        colvar_expval[i].is_derivative(); // constraints are not applied
        colvar_expval[i].reset();
      }
    }


    size_t t;
    std::string line;
    for (i = 0; i < colvars.size(); i++) {
       std::ifstream expfile (time_files[i]);
       if (expfile.is_open()) {
         for (t=0;t<colvar_centers[i].size();t++) {
            expfile >> colvar_times[i].vector1d_value[t];
            expfile >> colvar_expval[i].vector1d_value[t];
            getline (expfile,line);
         }
         expfile.close();
       } else {
         return cvm::error("Unable to open expFile", INPUT_ERROR);
       }
       colvar_centers[i]=colvar_expval[i];
    }

  }

  if (get_keyval(conf, "centers", colvar_centers, colvar_centers) && !read_exp_file) {
    for (i = 0; i < colvars.size(); i++) {
      if (cvm::debug()) {
        cvm::log("colvarbias_restraint: parsing initial centers, i = "+cvm::to_str(i)+".\n");
      }
      colvar_centers[i].apply_constraints();
    }
    null_centers = false;
  }

  if (null_centers && !read_exp_file) {
//    colvar_centers.clear();
    cvm::log("Warning: initial centers of the restraints set to zero (default choice for deerKernel cv).\n");
  }

  if (colvar_centers.size() != colvars.size()) {
    cvm::error("Error: number of centers does not match "
               "that of collective variables.\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  for (i = 0; i < colvars.size(); i++) {
     colvar_orig_centers[i]=colvar_centers[i]; // assign original centers
  }

  if (colvar_centers_errors.size() == 0) {
    colvar_centers_errors.resize(colvars.size());
    colvar_centers_errors.assign(colvars.size(), 0.0);
  }
  get_keyval(conf, "centersErrors", colvar_centers_errors, colvar_centers_errors);

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
  cvm::real eback;
  cvm::real lambdasum=0.0;
  cvm::real lambda2;
  cvm::real lambda2sum=0.0;
  if (opt_params) {
    cvm::real error_fact;
    colvar_cum_error=0;
    cvm::real cum_err;
    for (ii = 0; ii < numtypes[0]; ii++) {
      i=whichtypes[0][ii];
      cvm::real const unit_scale = 1.0/(colvar_widths_c[i] * colvar_widths_c[i]);
      error_fact = colvar_centers_errors[i]*colvar_centers_errors[i]/colvar_errors_scale;
      colvarvalue const deviation =
        0.5 * colvars[i]->dist2_lgrad(colvars[i]->value(), colvar_centers[i]);

      colvar_total_deviations[i] += weight * deviation * cvm::dt();
      bias_energy += kBT * unit_scale * colvar_total_deviations[i] * deviation;
      colvar_forces[i] = -1.0 * kBT * unit_scale * colvar_total_deviations[i];


      colvarvalue const error_drift = unit_scale * colvar_total_deviations[i]*error_fact;
      colvar_centers[i] = colvar_orig_centers[i] + error_drift;
      colvar_centers[i].apply_constraints();
      colvar_deviation[i]=error_drift/colvar_centers_errors[i];
      colvar_cum_error=colvar_cum_error+colvar_centers_errors[i]*colvar_centers_errors[i];
      for (t=0;t<colvar_centers[i].size();t++) {
         lambda2=-unit_scale * colvar_total_deviations[i].vector1d_value[t];
         lambda2=lambda2*lambda2;
         lambdasum+=colvar_centers_errors[i]*sqrt(lambda2);
         lambda2sum+=colvar_centers_errors[i]*colvar_centers_errors[i]*lambda2;;
      }
    }
    for (ii = 0; ii < numtypes[1]; ii++) {
      i=whichtypes[1][ii];
      cvm::real const unit_scale = 1.0/(colvar_widths_c[i] * colvar_widths_c[i]);
      colvarvalue const deviation =
        0.5 * colvars[i]->dist2_lgrad(colvars[i]->value(), colvar_centers[i]);

      colvar_total_deviations[i] += weight * deviation * cvm::dt();
      bias_energy += kBT * unit_scale * colvar_total_deviations[i] * deviation;
      colvar_forces[i] = -1.0 * kBT * unit_scale * colvar_total_deviations[i];
      cum_err=0.0;
      cvm::real coef_mdepth=0.0;
      cvm::real grad_mdepth=0.0;
      cvm::real coef_alpha=0.0;
      cvm::real grad_alpha=0.0;

      for (t=0;t<colvar_centers[i].size();t++) {
         eback=exp(-alpha_deer[ii]*sqrt(colvar_times[i].vector1d_value[t]*colvar_times[i].vector1d_value[t]));
         colvar_orig_centers[i].vector1d_value[t]=(colvar_expval[i].vector1d_value[t]-(1-mdepth_deer[ii])*eback)/(eback*mdepth_deer[ii]);
         error_fact =colvar_centers_errors[i]/(eback*mdepth_deer[ii]);
         cum_err=cum_err+error_fact*error_fact;
         error_fact = error_fact*error_fact/colvar_errors_scale;
         lambda2=-unit_scale * colvar_total_deviations[i].vector1d_value[t];
         cvm::real const drift=-lambda2*error_fact;
         colvar_centers[i].vector1d_value[t]=colvar_orig_centers[i].vector1d_value[t]+drift;
         colvar_deviation[i].vector1d_value[t]=drift*eback*mdepth_deer[ii]/colvar_centers_errors[i];
         //calculate derivatives of alpha_deer and mdepth_deer
         cvm::real const der_mdepth=(1-colvar_centers[i].vector1d_value[t]);
         coef_mdepth=coef_mdepth+(der_mdepth*der_mdepth);
         //coef_mdepth=coef_mdepth+(der_mdepth*(der_mdepth-(2*lambda2/unit_scale)));
         grad_mdepth=grad_mdepth+lambda2*der_mdepth;
         cvm::real const der_alpha=((1-mdepth_deer[ii])+mdepth_deer[ii]*colvar_centers[i].vector1d_value[t])*
                  sqrt(colvar_times[i].vector1d_value[t]*colvar_times[i].vector1d_value[t]);
         coef_alpha=coef_alpha+(der_alpha*der_alpha);
         //coef_alpha=coef_alpha+(der_alpha*(der_alpha+
         //          (mdepth_deer[ii]*lambda2*(sqrt(colvar_times[i].vector1d_value[t]*colvar_times[i].vector1d_value[t]))/unit_scale)));
         grad_alpha=grad_alpha+lambda2*der_alpha;
         lambda2=lambda2*lambda2;
         lambdasum+=(colvar_centers_errors[i]*sqrt(lambda2))/(mdepth_deer[ii]*eback);
         lambda2sum+=(colvar_centers_errors[i]*colvar_centers_errors[i]*lambda2)/(mdepth_deer[ii]*eback*mdepth_deer[ii]*eback);
      }
      cum_err=cum_err/colvar_centers[i].size();

      colvar_cum_error=colvar_cum_error+cum_err;
      colvar_centers[i].apply_constraints();
      // update alpha_deer and mdepth_deer
      coef_mdepth=sqrt(coef_mdepth*coef_mdepth);
      coef_alpha=sqrt(coef_alpha*coef_alpha);

      mdepth_deer[ii]=mdepth_deer[ii]-params_coupling_time * mdepth_deer[ii] * grad_mdepth * weight  * cvm::dt()/(unit_scale*coef_mdepth);
      alpha_deer[ii]=alpha_deer[ii]-params_coupling_time * mdepth_deer[ii] *  grad_alpha  * weight  * cvm::dt() /(unit_scale*coef_alpha);

    }
    colvar_cum_error=colvar_cum_error/double(colvars.size());
    if (use_norm_1) {
      for (i = 0; i < colvars.size(); i++) {
        cvm::real colvar_local_deviation=0.;
        colvar_deviation[i].set_absolute_value();
        colvarvalue const save_colvar_value = colvar_deviation[i];
        colvar_deviation[i].set_ones();
        colvar_local_deviation=colvar_deviation[i]*save_colvar_value;
        colvar_aver_deviation+=colvar_local_deviation;
      }
    } else {
      for (i = 0; i < colvars.size(); i++) {
        cvm::real colvar_local_deviation=0.;
        colvar_local_deviation=colvar_deviation[i]*colvar_deviation[i];
        colvar_aver_deviation+=colvar_local_deviation;
      }
    }
    // now update mdepth_deer and alpha_deer
  } else {
    for (i = 0; i < colvars.size(); i++) {
      cvm::real const unit_scale = 1.0/(colvar_widths_c[i] * colvar_widths_c[i]);
      cvm::real const error_fact = colvar_centers_errors[i]*colvar_centers_errors[i]/colvar_errors_scale;
      colvarvalue const deviation =
        0.5 * colvars[i]->dist2_lgrad(colvars[i]->value(), colvar_centers[i]);

      colvar_total_deviations[i] += weight * deviation * cvm::dt();
      bias_energy += kBT * unit_scale * colvar_total_deviations[i] * deviation;
      colvar_forces[i] = -1.0 * kBT * unit_scale * colvar_total_deviations[i];

      colvarvalue const error_drift = unit_scale * colvar_total_deviations[i]*error_fact;
      colvar_centers[i] = colvar_orig_centers[i] + error_drift;
      colvar_centers[i].apply_constraints();
      colvar_deviation[i]=error_drift/colvar_centers_errors[i];
      for (t=0;t<colvar_centers[i].size();t++) {
         lambda2=-unit_scale * colvar_total_deviations[i].vector1d_value[t];
         lambda2=lambda2*lambda2;
         lambdasum+=colvar_centers_errors[i]*sqrt(lambda2);
         lambda2sum+=colvar_centers_errors[i]*colvar_centers_errors[i]*lambda2;
      }

    }
    if (use_norm_1) {
      for (i = 0; i < colvars.size(); i++) {
        cvm::real colvar_local_deviation=0.;
        colvar_deviation[i].set_absolute_value();
        colvarvalue const save_colvar_value = colvar_deviation[i];
        colvar_deviation[i].set_ones();
        colvar_local_deviation=colvar_deviation[i]*save_colvar_value;
        colvar_aver_deviation+=colvar_local_deviation;
      }
    } else {
      for (i = 0; i < colvars.size(); i++) {
        cvm::real colvar_local_deviation=0.;
        colvar_local_deviation=colvar_deviation[i]*colvar_deviation[i];
        colvar_aver_deviation+=colvar_local_deviation;
      }
    }
  }


  colvar_aver_deviation=colvar_aver_deviation/double(colvar_size_tot);

  lambdasum=lambdasum/double(colvar_size_tot);
  lambda2sum=sqrt(lambda2sum/double(colvar_size_tot));

  if (!use_norm_1) {
    colvar_aver_deviation=sqrt(colvar_aver_deviation);
  }

  if (fix_chi_square_one) {

//    printf("FIX ONE: %i %30.20f %30.20f %30.20f %30.20f %30.20f \n",cvm::step_absolute(),colvar_aver_deviation,colvar_cum_error,colvar_cum_uscale,colvar_errors_scale,double(colvar_size_tot));
    colvar_errors_scale=lambda2sum;
    if (use_norm_1) colvar_errors_scale=lambdasum;
    if (colvar_errors_scale<1.0/colvar_maxerror_scale) colvar_errors_scale=1.0/colvar_maxerror_scale;
    if (colvar_errors_scale>colvar_maxerror_scale) colvar_errors_scale=colvar_maxerror_scale;

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

    if (opt_params) {
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
         << sqrt(colvar_cum_error/colvar_errors_scale)
         << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_size_tot
         << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_cum_error
         << "\n";
      error_code |= (os.good() ? COLVARS_OK : FILE_ERROR);
    }

    if (rad_param_os) {
      std::ostream &os = *rad_param_os;
      os.setf(std::ios::scientific, std::ios::floatfield);
      os << std::setw(cvm::it_width) << cvm::step_absolute()
         << cvm::step_absolute()
         << " ";
      for (int ii = 0; ii < numtypes[1]; ii++) {
        os << " "
           << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
           << mdepth_deer[ii]
           << " "
           << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
           << alpha_deer[ii];
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
  for (size_t i = 0; i < colvars.size(); i++) {
    os << colvar_total_deviations[i];
  }
  if (fix_chi_square_one) {
    os << "    error_scale " << colvar_errors_scale << "\n";
  }

  if (opt_params) {
    if (numtypes[1]>0) {
      os << "    deerMdepth ";
      for (size_t i = 0; i < numtypes[1]; i++) {
         os << mdepth_deer[i];
         os << " ";
      }
      os << "\n";
      os << "    deerBackAlpha ";
      for (size_t i = 0; i < numtypes[1]; i++) {
         os << alpha_deer[i];
         os << " ";
      }
      os << "\n";
    }
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

  if (opt_params) {
    // read last value of the parameters
    if (numtypes[1]>0) {
      if (!get_keyval(state_conf, "deerMdepth", mdepth_deer))
        cvm::error("Error: missing deerMdepth from the restart.\n");
      if (!get_keyval(state_conf, "deerBackAlpha", alpha_deer))
        cvm::error("Error: missing deerBackAlpha from the restart.\n");
    }
  }

  return error_code;
}
