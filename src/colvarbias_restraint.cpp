// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias_restraint.h"



colvarbias_restraint::colvarbias_restraint(char const *key)
  : colvarbias(key)
{
}


int colvarbias_restraint::init(std::string const &conf)
{
  colvarbias::init(conf);
  enable(f_cvb_apply_force);

  if (cvm::debug())
    cvm::log("Initializing a new restraint bias.\n");

  return COLVARS_OK;
}


int colvarbias_restraint::update()
{
  bias_energy = 0.0;

  if (cvm::debug())
    cvm::log("Updating the restraint bias \""+this->name+"\".\n");

  // Force and energy calculation
  for (size_t i = 0; i < colvars.size(); i++) {
    colvar_forces[i].type(colvars[i]->value());
    colvar_forces[i].is_derivative();
    colvar_forces[i] = restraint_force(i);
    bias_energy += restraint_potential(i);
  }

  if (cvm::debug())
    cvm::log("Done updating the restraint bias \""+this->name+"\".\n");

  if (cvm::debug())
    cvm::log("Current forces for the restraint bias \""+
             this->name+"\": "+cvm::to_str(colvar_forces)+".\n");

  return COLVARS_OK;
}


colvarbias_restraint::~colvarbias_restraint()
{
  if (cvm::n_rest_biases > 0)
    cvm::n_rest_biases -= 1;
}


std::string const colvarbias_restraint::get_state_params() const
{
  return colvarbias::get_state_params();
}


int colvarbias_restraint::set_state_params(std::string const &conf)
{
  return colvarbias::set_state_params(conf);
}


std::ostream & colvarbias_restraint::write_traj_label(std::ostream &os)
{
  return colvarbias::write_traj_label(os);
}


std::ostream & colvarbias_restraint::write_traj(std::ostream &os)
{
  return colvarbias::write_traj(os);
}



colvarbias_restraint_centers::colvarbias_restraint_centers(char const *key)
  : colvarbias(key), colvarbias_restraint(key)
{
}


int colvarbias_restraint_centers::init(std::string const &conf)
{
  size_t i;

  bool null_centers = (colvar_centers.size() == 0);
  if (null_centers) {
    // try to initialize the restraint centers for the first time
    colvar_centers.resize(colvars.size());
    colvar_centers_raw.resize(colvars.size());
    for (i = 0; i < colvars.size(); i++) {
      colvar_centers[i].type(colvars[i]->value());
      colvar_centers[i].reset();
      colvar_centers_raw[i].type(colvars[i]->value());
      colvar_centers_raw[i].reset();
    }
  }

  if (get_keyval(conf, "centers", colvar_centers, colvar_centers)) {
    for (i = 0; i < colvars.size(); i++) {
      if (cvm::debug()) {
        cvm::log("colvarbias_restraint: parsing initial centers, i = "+cvm::to_str(i)+".\n");
      }
      colvar_centers_raw[i] = colvar_centers[i];
      colvar_centers[i].apply_constraints();
    }
    null_centers = false;
  }

  if (null_centers) {
    colvar_centers.clear();
    cvm::error("Error: must define the initial centers of the restraints.\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  if (colvar_centers.size() != colvars.size()) {
    cvm::error("Error: number of centers does not match "
               "that of collective variables.\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  return COLVARS_OK;
}



colvarbias_restraint_k::colvarbias_restraint_k(char const *key)
  : colvarbias(key), colvarbias_restraint(key)
{
  force_k = 1.0;
}


int colvarbias_restraint_k::init(std::string const &conf)
{
  get_keyval(conf, "forceConstant", force_k, force_k);
  return COLVARS_OK;
}



colvarbias_restraint_moving::colvarbias_restraint_moving(char const *key)
{
  target_nstages = 0;
  target_nsteps = 0;
  stage = 0;
  b_chg_centers = false;
  b_chg_force_k = false;
}


int colvarbias_restraint_moving::init(std::string const &conf)
{
  if (b_chg_centers && b_chg_force_k) {
    cvm::error("Error: cannot specify both targetCenters and targetForceConstant.\n",
               INPUT_ERROR);
    return INPUT_ERROR;
  }

  if (b_chg_centers || b_chg_force_k) {

    get_keyval(conf, "targetNumSteps", target_nsteps, target_nsteps);
    if (!target_nsteps) {
      cvm::error("Error: targetNumSteps must be non-zero.\n", INPUT_ERROR);
      return cvm::get_error();
    }

    if (get_keyval(conf, "targetNumStages", target_nstages, target_nstages) &&
        lambda_schedule.size()) {
      cvm::error("Error: targetNumStages and lambdaSchedule are incompatible.\n", INPUT_ERROR);
      return cvm::get_error();
    }
  }

  return COLVARS_OK;
}


std::string const colvarbias_restraint_moving::get_state_params() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  if (b_chg_centers || b_chg_force_k) {
    // TODO move this
    if (target_nstages) {
      os << "stage " << std::setw(cvm::it_width)
         << stage << "\n";
    }
  }
  return os.str();
}


int colvarbias_restraint_moving::set_state_params(std::string const &conf)
{
  if (b_chg_centers || b_chg_force_k) {
    if (target_nstages) {
      //    cvm::log ("Reading current stage from the restart.\n");
      if (!get_keyval(conf, "stage", stage))
        cvm::error("Error: current stage is missing from the restart.\n");
    }
  }
  return COLVARS_OK;
}



colvarbias_restraint_centers_moving::colvarbias_restraint_centers_moving(char const *key)
  : colvarbias(key),
    colvarbias_restraint(key),
    colvarbias_restraint_moving(key),
    colvarbias_restraint_centers(key)
{
  b_chg_centers = false;
  b_output_centers = false;
  b_output_acc_work = false;
  acc_work = 0.0;
}


int colvarbias_restraint_centers_moving::init(std::string const &conf)
{
  colvarbias_restraint_centers::init(conf);

  if (cvm::debug()) {
    cvm::log("colvarbias_restraint: parsing target centers.\n");
  }

  size_t i;
  if (get_keyval(conf, "targetCenters", target_centers, colvar_centers)) {
    if (colvar_centers.size() != colvars.size()) {
      cvm::error("Error: number of target centers does not match "
                 "that of collective variables.\n");
    }
    b_chg_centers = true;
    for (i = 0; i < target_centers.size(); i++) {
      target_centers[i].apply_constraints();
    }
  }

  if (b_chg_centers) {
    // parse moving restraint options
    colvarbias_restraint_moving::init(conf);
  } else {
    target_centers.clear();
    return COLVARS_OK;
  }

  get_keyval(conf, "outputCenters", b_output_centers, b_output_centers);
  get_keyval(conf, "outputAccumulatedWork", b_output_acc_work, b_output_acc_work);

  return COLVARS_OK;
}


int colvarbias_restraint_centers_moving::update()
{
  if (b_chg_centers) {

    if (cvm::debug()) {
      cvm::log("Updating centers for the restraint bias \""+
               this->name+"\": "+cvm::to_str(colvar_centers)+".\n");
    }

    if (!centers_incr.size()) {
      // if this is the first calculation, calculate the advancement
      // at each simulation step (or stage, if applicable)
      // (take current stage into account: it can be non-zero
      //  if we are restarting a staged calculation)
      centers_incr.resize(colvars.size());
      for (size_t i = 0; i < colvars.size(); i++) {
        centers_incr[i].type(colvars[i]->value());
        centers_incr[i] = (target_centers[i] - colvar_centers_raw[i]) /
          cvm::real( target_nstages ? (target_nstages - stage) :
                     (target_nsteps - cvm::step_absolute()));
      }
      if (cvm::debug()) {
        cvm::log("Center increment for the restraint bias \""+
                 this->name+"\": "+cvm::to_str(centers_incr)+" at stage "+cvm::to_str(stage)+ ".\n");
      }
    }

    if (target_nstages) {
      if ((cvm::step_relative() > 0)
          && (cvm::step_absolute() % target_nsteps) == 0
          && stage < target_nstages) {

        for (size_t i = 0; i < colvars.size(); i++) {
          colvar_centers_raw[i] += centers_incr[i];
          colvar_centers[i] = colvar_centers_raw[i];
          colvars[i]->wrap(colvar_centers[i]);
          colvar_centers[i].apply_constraints();
        }
        stage++;
        cvm::log("Moving restraint \"" + this->name +
                 "\" stage " + cvm::to_str(stage) +
                 " : setting centers to " + cvm::to_str(colvar_centers) +
                 " at step " +  cvm::to_str(cvm::step_absolute()));
      }
    } else if ((cvm::step_relative() > 0) && (cvm::step_absolute() <= target_nsteps)) {
      // move the restraint centers in the direction of the targets
      // (slow growth)
      for (size_t i = 0; i < colvars.size(); i++) {
        colvar_centers_raw[i] += centers_incr[i];
        colvar_centers[i] = colvar_centers_raw[i];
        colvars[i]->wrap(colvar_centers[i]);
        colvar_centers[i].apply_constraints();
      }
    }

    if (cvm::debug()) {
      cvm::log("New centers for the restraint bias \""+
               this->name+"\": "+cvm::to_str(colvar_centers)+".\n");
    }
  }

  return COLVARS_OK;
}


int colvarbias_restraint_centers_moving::update_acc_work()
{
  if (b_output_acc_work) {
    if ((cvm::step_relative() > 0) || (cvm::step_absolute() == 0)) {
      for (size_t i = 0; i < colvars.size(); i++) {
        // project forces on the calculated increments at this step
        acc_work += colvar_forces[i] * centers_incr[i];
      }
    }
  }
  return COLVARS_OK;
}


std::string const colvarbias_restraint_centers_moving::get_state_params() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);

  if (b_chg_centers) {
    size_t i;
    os << "centers ";
    for (i = 0; i < colvars.size(); i++) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_centers[i];
    }
    os << "\n";
    os << "centers_raw ";
    for (i = 0; i < colvars.size(); i++) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_centers_raw[i];
    }
    os << "\n";

    if (b_output_acc_work) {
      os << "accumulatedWork "
         << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width)
         << acc_work << "\n";
    }
  }

  return colvarbias_restraint_moving::get_state_params() + os.str();
}


int colvarbias_restraint_centers_moving::set_state_params(std::string const &conf)
{
  colvarbias_restraint::set_state_params(conf);

  if (b_chg_centers) {
    //    cvm::log ("Reading the updated restraint centers from the restart.\n");
    if (!get_keyval(conf, "centers", colvar_centers))
      cvm::error("Error: restraint centers are missing from the restart.\n");
    if (!get_keyval(conf, "centers_raw", colvar_centers_raw))
      cvm::error("Error: \"raw\" restraint centers are missing from the restart.\n");
    if (b_output_acc_work) {
      if (!get_keyval(conf, "accumulatedWork", acc_work))
        cvm::error("Error: accumulatedWork is missing from the restart.\n");
    }
  }

  return COLVARS_OK;
}


std::ostream & colvarbias_restraint_centers_moving::write_traj_label(std::ostream &os)
{
  if (b_output_centers) {
    for (size_t i = 0; i < colvars.size(); i++) {
      size_t const this_cv_width = (colvars[i]->value()).output_width(cvm::cv_width);
      os << " x0_"
         << cvm::wrap_string(colvars[i]->name, this_cv_width-3);
    }
  }

  if (b_output_acc_work) {
    os << " W_"
       << cvm::wrap_string(this->name, cvm::en_width-2);
  }

  return os;
}


std::ostream & colvarbias_restraint_centers_moving::write_traj(std::ostream &os)
{
  if (b_output_centers) {
    for (size_t i = 0; i < colvars.size(); i++) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_centers[i];
    }
  }

  if (b_output_acc_work) {
    os << " "
       << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width)
       << acc_work;
  }

  return os;
}



colvarbias_restraint_k_moving::colvarbias_restraint_k_moving(char const *key)
  : colvarbias(key),
    colvarbias_restraint(key),
    colvarbias_restraint_k(key),
    colvarbias_restraint_moving(key)
{
  b_chg_force_k = false;
  target_equil_steps = 0;
  target_force_k = 0.0;
  starting_force_k = 0.0;
  force_k_exp = 1.0;
  restraint_FE = 0.0;
}


int colvarbias_restraint_k_moving::init(std::string const &conf)
{
  colvarbias_restraint_k::init(conf);

  if (get_keyval(conf, "targetForceConstant", target_force_k, target_force_k)) {
    starting_force_k = force_k;
    b_chg_force_k = true;
  }

  if (b_chg_force_k) {
    // parse moving restraint options
    colvarbias_restraint_moving::init(conf);
  } else {
    return COLVARS_OK;
  }

  get_keyval(conf, "targetEquilSteps", target_equil_steps, target_equil_steps);

  get_keyval(conf, "lambdaSchedule", lambda_schedule, lambda_schedule);
  if (lambda_schedule.size()) {
    // There is one more lambda-point than stages
    target_nstages = lambda_schedule.size() - 1;
  }

  if (get_keyval(conf, "targetForceExponent", force_k_exp, force_k_exp)) {
    if (! b_chg_force_k)
      cvm::log("Warning: not changing force constant: targetForceExponent will be ignored\n");
  }
  if (force_k_exp < 1.0) {
    cvm::log("Warning: for all practical purposes, targetForceExponent should be 1.0 or greater.\n");
  }

  return COLVARS_OK;
}


int colvarbias_restraint_k_moving::update()
{
  if (b_chg_force_k) {

    cvm::real lambda;

    if (target_nstages) {

      if (cvm::step_absolute() == 0) {
        // Setup first stage of staged variable force constant calculation
        if (lambda_schedule.size()) {
          lambda = lambda_schedule[0];
        } else {
          lambda = 0.0;
        }
        force_k = starting_force_k + (target_force_k - starting_force_k)
          * std::pow(lambda, force_k_exp);
        cvm::log("Restraint " + this->name + ", stage " +
                 cvm::to_str(stage) + " : lambda = " + cvm::to_str(lambda));
        cvm::log("Setting force constant to " + cvm::to_str(force_k));
      }

      // TI calculation: estimate free energy derivative
      // need current lambda
      if (lambda_schedule.size()) {
        lambda = lambda_schedule[stage];
      } else {
        lambda = cvm::real(stage) / cvm::real(target_nstages);
      }

      if (target_equil_steps == 0 || cvm::step_absolute() % target_nsteps >= target_equil_steps) {
        // Start averaging after equilibration period, if requested

        // Square distance normalized by square colvar width
        cvm::real dist_sq = 0.0;
        for (size_t i = 0; i < colvars.size(); i++) {
          dist_sq += d_restraint_potential_dk(i);
        }

        restraint_FE += 0.5 * force_k_exp * std::pow(lambda, force_k_exp - 1.0)
          * (target_force_k - starting_force_k) * dist_sq;
      }

      // Finish current stage...
      if (cvm::step_absolute() % target_nsteps == 0 &&
          cvm::step_absolute() > 0) {

        cvm::log("Lambda= " + cvm::to_str(lambda) + " dA/dLambda= "
                 + cvm::to_str(restraint_FE / cvm::real(target_nsteps - target_equil_steps)));

        //  ...and move on to the next one
        if (stage < target_nstages) {

          restraint_FE = 0.0;
          stage++;
          if (lambda_schedule.size()) {
            lambda = lambda_schedule[stage];
          } else {
            lambda = cvm::real(stage) / cvm::real(target_nstages);
          }
          force_k = starting_force_k + (target_force_k - starting_force_k)
            * std::pow(lambda, force_k_exp);
          cvm::log("Restraint " + this->name + ", stage " +
                   cvm::to_str(stage) + " : lambda = " + cvm::to_str(lambda));
          cvm::log("Setting force constant to " + cvm::to_str(force_k));
        }
      }

    } else if (cvm::step_absolute() <= target_nsteps) {

      // update force constant (slow growth)
      lambda = cvm::real(cvm::step_absolute()) / cvm::real(target_nsteps);
      force_k = starting_force_k + (target_force_k - starting_force_k)
        * std::pow(lambda, force_k_exp);
    }
  }

  return COLVARS_OK;
}


std::string const colvarbias_restraint_k_moving::get_state_params() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  if (b_chg_force_k) {
    os << "forceConstant "
       << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << force_k << "\n";
  }
  return colvarbias_restraint_moving::get_state_params() + os.str();
}


int colvarbias_restraint_k_moving::set_state_params(std::string const &conf)
{
  colvarbias_restraint::set_state_params(conf);

  if (b_chg_force_k) {
    //    cvm::log ("Reading the updated force constant from the restart.\n");
    if (!get_keyval(conf, "forceConstant", force_k, force_k))
      cvm::error("Error: force constant is missing from the restart.\n");
  }

  return COLVARS_OK;
}


std::ostream & colvarbias_restraint_k_moving::write_traj_label(std::ostream &os)
{
  return os;
}


std::ostream & colvarbias_restraint_k_moving::write_traj(std::ostream &os)
{
  return os;
}


// // TODO remove these two
// void colvarbias_restraint::change_configuration(std::string const &conf)
// {
//   get_keyval(conf, "forceConstant", force_k, force_k);
//   if (get_keyval(conf, "centers", colvar_centers, colvar_centers)) {
//     for (size_t i = 0; i < colvars.size(); i++) {
//       colvar_centers[i].type(colvars[i]->value());
//       colvar_centers[i].apply_constraints();
//       colvar_centers_raw[i].type(colvars[i]->value());
//       colvar_centers_raw[i] = colvar_centers[i];
//     }
//   }
// }


// cvm::real colvarbias_restraint::energy_difference(std::string const &conf)
// {
//   std::vector<colvarvalue> alt_colvar_centers;
//   cvm::real alt_force_k;
//   cvm::real alt_bias_energy = 0.0;

//   get_keyval(conf, "forceConstant", alt_force_k, force_k);

//   alt_colvar_centers.resize(colvars.size());
//   size_t i;
//   for (i = 0; i < colvars.size(); i++) {
//     alt_colvar_centers[i].type(colvars[i]->value());
//   }
//   if (get_keyval(conf, "centers", alt_colvar_centers, colvar_centers)) {
//     for (i = 0; i < colvars.size(); i++) {
//       alt_colvar_centers[i].apply_constraints();
//     }
//   }

//   for (i = 0; i < colvars.size(); i++) {
//     alt_bias_energy += restraint_potential(alt_force_k,
// 					   colvars[i],
// 					   alt_colvar_centers[i]);
//   }

//   return alt_bias_energy - bias_energy;
// }





// redefined due to legacy state file keyword "harmonic"
std::istream & colvarbias_restraint::read_state(std::istream &is)
{
  size_t const start_pos = is.tellg();

  std::string key, brace, conf;
  if ( !(is >> key)   || !(key == "restraint" || key == "harmonic") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block("configuration", conf)) ||
       (set_state_params(conf) != COLVARS_OK) ) {
    cvm::error("Error: in reading state configuration for \""+bias_type+"\" bias \""+
               this->name+"\" at position "+
               cvm::to_str(is.tellg())+" in stream.\n", INPUT_ERROR);
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  if (!read_state_data(is)) {
    cvm::error("Error: in reading state data for \""+bias_type+"\" bias \""+
               this->name+"\" at position "+
               cvm::to_str(is.tellg())+" in stream.\n", INPUT_ERROR);
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
  }

  is >> brace;
  if (brace != "}") {
    cvm::log("brace = "+brace+"\n");
    cvm::error("Error: corrupt restart information for \""+bias_type+"\" bias \""+
               this->name+"\": no matching brace at position "+
               cvm::to_str(is.tellg())+" in stream.\n");
    is.setstate(std::ios::failbit);
  }

  return is;
}


std::ostream & colvarbias_restraint::write_state(std::ostream &os)
{
  os.setf(std::ios::scientific, std::ios::floatfield);
  os << "restraint {\n"
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



colvarbias_restraint_harmonic::colvarbias_restraint_harmonic(char const *key)
  : colvarbias(key),
    colvarbias_restraint(key),
    colvarbias_restraint_centers(key),
    colvarbias_restraint_k(key),
    colvarbias_restraint_moving(key),
    colvarbias_restraint_centers_moving(key),
    colvarbias_restraint_k_moving(key)
{
}


int colvarbias_restraint_harmonic::init(std::string const &conf)
{
  colvarbias_restraint::init(conf);
  colvarbias_restraint_moving::init(conf);
  colvarbias_restraint_centers_moving::init(conf);
  colvarbias_restraint_k_moving::init(conf);

  for (size_t i = 0; i < colvars.size(); i++) {
    if (colvars[i]->width != 1.0)
      cvm::log("The force constant for colvar \""+colvars[i]->name+
               "\" will be rescaled to "+
               cvm::to_str(force_k / (colvars[i]->width * colvars[i]->width))+
               " according to the specified width.\n");
  }

  return COLVARS_OK;
}


int colvarbias_restraint_harmonic::update()
{
  // update parameters (centers or force constant)
  colvarbias_restraint_centers_moving::update();
  colvarbias_restraint_k_moving::update();

  // update restraint energy and forces
  colvarbias_restraint::update();

  // update accumulated work using the current forces
  colvarbias_restraint_centers_moving::update_acc_work();

  return COLVARS_OK;
}


cvm::real colvarbias_restraint_harmonic::restraint_potential(size_t i) const
{
  return 0.5 * force_k / (colvars[i]->width * colvars[i]->width) *
    colvars[i]->dist2(colvars[i]->value(), colvar_centers[i]);
}


colvarvalue const colvarbias_restraint_harmonic::restraint_force(size_t i) const
{
  return -0.5 * force_k / (colvars[i]->width * colvars[i]->width) *
    colvars[i]->dist2_lgrad(colvars[i]->value(), colvar_centers[i]);
}


cvm::real colvarbias_restraint_harmonic::d_restraint_potential_dk(size_t i) const
{
  return 0.5 / (colvars[i]->width * colvars[i]->width) *
    colvars[i]->dist2(colvars[i]->value(), colvar_centers[i]);
}


std::string const colvarbias_restraint_harmonic::get_state_params() const
{
  return colvarbias_restraint::get_state_params() +
    colvarbias_restraint_centers_moving::get_state_params() +
    colvarbias_restraint_k_moving::get_state_params();
}


int colvarbias_restraint_harmonic::set_state_params(std::string const &conf)
{
  int error_code = COLVARS_OK;
  cvm::combine_errors(error_code, colvarbias_restraint::set_state_params(conf));
  cvm::combine_errors(error_code, colvarbias_restraint_centers_moving::set_state_params(conf));
  cvm::combine_errors(error_code, colvarbias_restraint_k_moving::set_state_params(conf));
  return error_code;
}


std::ostream & colvarbias_restraint_harmonic::write_traj_label(std::ostream &os)
{
  colvarbias_restraint::write_traj_label(os);
  colvarbias_restraint_centers_moving::write_traj_label(os);
  colvarbias_restraint_k_moving::write_traj_label(os);
  return os;
}


std::ostream & colvarbias_restraint_harmonic::write_traj(std::ostream &os)
{
  colvarbias_restraint::write_traj(os);
  colvarbias_restraint_centers_moving::write_traj(os);
  colvarbias_restraint_k_moving::write_traj(os);
  return os;
}



colvarbias_restraint_linear::colvarbias_restraint_linear(char const *key)
  : colvarbias(key),
    colvarbias_restraint(key),
    colvarbias_restraint_centers(key),
    colvarbias_restraint_k(key),
    colvarbias_restraint_moving(key),
    colvarbias_restraint_centers_moving(key),
    colvarbias_restraint_k_moving(key)
{
}


int colvarbias_restraint_linear::init(std::string const &conf)
{
  colvarbias_restraint::init(conf);
  colvarbias_restraint_moving::init(conf);
  colvarbias_restraint_centers_moving::init(conf);
  colvarbias_restraint_k_moving::init(conf);

  for (size_t i = 0; i < colvars.size(); i++) {
    if (colvars[i]->is_enabled(f_cv_periodic)) {
      cvm::error("Error: linear biases cannot be applied to periodic variables.\n",
                 INPUT_ERROR);
      return INPUT_ERROR;
    }
    if (colvars[i]->width != 1.0)
      cvm::log("The force constant for colvar \""+colvars[i]->name+
               "\" will be rescaled to "+
               cvm::to_str(force_k / colvars[i]->width)+
               " according to the specified width.\n");
  }

  return COLVARS_OK;
}


int colvarbias_restraint_linear::update()
{
  // update parameters (centers or force constant)
  colvarbias_restraint_centers_moving::update();
  colvarbias_restraint_k_moving::update();

  // update restraint energy and forces
  colvarbias_restraint::update();

  // update accumulated work using the current forces
  colvarbias_restraint_centers_moving::update_acc_work();

  return COLVARS_OK;
}


cvm::real colvarbias_restraint_linear::restraint_potential(size_t i) const
{
  return force_k / colvars[i]->width * (colvars[i]->value() - colvar_centers[i]);
}


colvarvalue const colvarbias_restraint_linear::restraint_force(size_t i) const
{
  return -1.0 * force_k / colvars[i]->width;
}


cvm::real colvarbias_restraint_linear::d_restraint_potential_dk(size_t i) const
{
  return 1.0 / colvars[i]->width * (colvars[i]->value() - colvar_centers[i]);
}


std::string const colvarbias_restraint_linear::get_state_params() const
{
  return colvarbias_restraint::get_state_params() +
    colvarbias_restraint_centers_moving::get_state_params() +
    colvarbias_restraint_k_moving::get_state_params();
}


int colvarbias_restraint_linear::set_state_params(std::string const &conf)
{
  int error_code = COLVARS_OK;
  cvm::combine_errors(error_code, colvarbias_restraint::set_state_params(conf));
  cvm::combine_errors(error_code, colvarbias_restraint_centers_moving::set_state_params(conf));
  cvm::combine_errors(error_code, colvarbias_restraint_k_moving::set_state_params(conf));
  return error_code;
}


std::ostream & colvarbias_restraint_linear::write_traj_label(std::ostream &os)
{
  colvarbias_restraint::write_traj_label(os);
  colvarbias_restraint_centers_moving::write_traj_label(os);
  colvarbias_restraint_k_moving::write_traj_label(os);
  return os;
}


std::ostream & colvarbias_restraint_linear::write_traj(std::ostream &os)
{
  colvarbias_restraint::write_traj(os);
  colvarbias_restraint_centers_moving::write_traj(os);
  colvarbias_restraint_k_moving::write_traj(os);
  return os;
}
