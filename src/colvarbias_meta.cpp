// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <fstream>
#include <iomanip>
#include <algorithm>

// Define function to get the absolute path of a replica file
#if defined(_WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define GETCWD(BUF, SIZE) ::_getcwd(BUF, SIZE)
#define PATHSEP "\\"
#else
#include <unistd.h>
#define GETCWD(BUF, SIZE) ::getcwd(BUF, SIZE)
#define PATHSEP "/"
#endif

#ifdef __cpp_lib_filesystem
// When std::filesystem is available, use it
#include <filesystem>
#undef GETCWD
#define GETCWD(BUF, SIZE) (std::filesystem::current_path().string().c_str())
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias_meta.h"
#include "colvars_memstream.h"


colvarbias_meta::colvarbias_meta(char const *key)
  : colvarbias(key), colvarbias_ti(key)
{
  new_hills_begin = hills.end();

  hill_weight = 0.0;
  hill_width = 0.0;

  new_hill_freq = 1000;

  use_grids = true;
  grids_freq = 0;
  rebin_grids = false;
  hills_energy = NULL;
  hills_energy_gradients = NULL;

  dump_fes = true;
  keep_hills = false;
  restart_keep_hills = false;
  dump_fes_save = false;
  dump_replica_fes = false;

  b_hills_traj = false;

  ebmeta_equil_steps = 0L;

  replica_id.clear();
}


int colvarbias_meta::init(std::string const &conf)
{
  int error_code = COLVARS_OK;
  size_t i = 0;

  error_code |= colvarbias::init(conf);
  error_code |= colvarbias_ti::init(conf);

  cvm::main()->cite_feature("Metadynamics colvar bias implementation");

  enable(f_cvb_calc_pmf);

  get_keyval(conf, "hillWeight", hill_weight, hill_weight);
  if (hill_weight > 0.0) {
    enable(f_cvb_apply_force);
  } else {
    cvm::error("Error: hillWeight must be provided, and a positive number.\n", COLVARS_INPUT_ERROR);
  }

  get_keyval(conf, "newHillFrequency", new_hill_freq, new_hill_freq);
  if (new_hill_freq > 0) {
    enable(f_cvb_history_dependent);
    if (grids_freq == 0) {
      grids_freq = new_hill_freq;
    }
  }

  get_keyval(conf, "gaussianSigmas", colvar_sigmas, colvar_sigmas);

  get_keyval(conf, "hillWidth", hill_width, hill_width);

  if ((colvar_sigmas.size() > 0) && (hill_width > 0.0)) {
    error_code |= cvm::error("Error: hillWidth and gaussianSigmas are "
                             "mutually exclusive.", COLVARS_INPUT_ERROR);
  }

  if (hill_width > 0.0) {
    colvar_sigmas.resize(num_variables());
    // Print the calculated sigma parameters
    cvm::log("Half-widths of the Gaussian hills (sigma's):\n");
    for (i = 0; i < num_variables(); i++) {
      colvar_sigmas[i] = variables(i)->width * hill_width / 2.0;
      cvm::log(variables(i)->name+std::string(": ")+
               cvm::to_str(colvar_sigmas[i]));
    }
  }

  if (colvar_sigmas.size() == 0) {
    error_code |= cvm::error("Error: positive values are required for "
                             "either hillWidth or gaussianSigmas.",
                             COLVARS_INPUT_ERROR);
  }

  {
    bool b_replicas = false;
    get_keyval(conf, "multipleReplicas", b_replicas, false);
    if (b_replicas) {
      cvm::main()->cite_feature("Multiple-walker metadynamics colvar bias implementation");
  comm = multiple_replicas;
    } else {
      comm = single_replica;
    }
  }

  get_keyval(conf, "useGrids", use_grids, use_grids);

  if (use_grids) {

    for (i = 0; i < num_variables(); i++) {
      if (2.0*colvar_sigmas[i] < variables(i)->width) {
        cvm::log("Warning: gaussianSigmas is too narrow for the grid "
                 "spacing along "+variables(i)->name+".");
      }
    }

    get_keyval(conf, "gridsUpdateFrequency", grids_freq, grids_freq);
    get_keyval(conf, "rebinGrids", rebin_grids, rebin_grids);

    expand_grids = false;
    for (i = 0; i < num_variables(); i++) {
      variables(i)->enable(f_cv_grid); // Could be a child dependency of a f_cvb_use_grids feature
      if (variables(i)->expand_boundaries) {
        expand_grids = true;
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                 ": Will expand grids when the colvar \""+
                 variables(i)->name+"\" approaches its boundaries.\n");
      }
    }

    get_keyval(conf, "writeFreeEnergyFile", dump_fes, dump_fes);

    get_keyval(conf, "keepHills", keep_hills, keep_hills);
    get_keyval(conf, "keepFreeEnergyFiles", dump_fes_save, dump_fes_save);

    if (hills_energy == NULL) {
      hills_energy           = new colvar_grid_scalar(colvars);
      hills_energy_gradients = new colvar_grid_gradient(colvars);
    }

  } else {

    dump_fes = false;
  }

  get_keyval(conf, "writeHillsTrajectory", b_hills_traj, b_hills_traj);

  error_code |= init_replicas_params(conf);
  error_code |= init_well_tempered_params(conf);
  error_code |= init_reflection_params(conf);
  error_code |= init_interval_params(conf);
  error_code |= init_ebmeta_params(conf);

  if (cvm::debug())
    cvm::log("Done initializing the metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+".\n");

  return error_code;
}


int colvarbias_meta::init_replicas_params(std::string const &conf)
{
  colvarproxy *proxy = cvm::main()->proxy;

  // in all cases, the first replica is this bias itself
  if (replicas.size() == 0) {
    replicas.push_back(this);
  }

  if (comm != single_replica) {

    if (!get_keyval(conf, "writePartialFreeEnergyFile",
                    dump_replica_fes, dump_replica_fes)) {
      get_keyval(conf, "dumpPartialFreeEnergyFile", dump_replica_fes,
                 dump_replica_fes, colvarparse::parse_silent);
    }

    if (dump_replica_fes && (! dump_fes)) {
      dump_fes = true;
      cvm::log("Enabling \"writeFreeEnergyFile\".\n");
    }

    get_keyval(conf, "replicaID", replica_id, replica_id);
    if (!replica_id.size()) {
      if (proxy->replica_enabled() == COLVARS_OK) {
        // Obtain replicaID from the communicator
        replica_id = cvm::to_str(proxy->replica_index());
        cvm::log("Setting replicaID from communication layer: replicaID = "+
                 replica_id+".\n");
      } else {
        return cvm::error("Error: using more than one replica, but replicaID "
                          "could not be obtained.\n", COLVARS_INPUT_ERROR);
      }
    }

    get_keyval(conf, "replicasRegistry", replicas_registry_file,
               replicas_registry_file);
    if (!replicas_registry_file.size()) {
      return cvm::error("Error: the name of the \"replicasRegistry\" file "
                        "must be provided.\n", COLVARS_INPUT_ERROR);
    }

    get_keyval(conf, "replicaUpdateFrequency",
               replica_update_freq, replica_update_freq);
    if (replica_update_freq == 0) {
      return cvm::error("Error: replicaUpdateFrequency must be positive.\n",
                        COLVARS_INPUT_ERROR);
    }

    if (expand_grids) {
      return cvm::error("Error: expandBoundaries is not supported when "
                        "using more than one replicas; please allocate "
                        "wide enough boundaries for each colvar"
                        "ahead of time.\n", COLVARS_INPUT_ERROR);
    }

    if (keep_hills) {
      return cvm::error("Error: multipleReplicas and keepHills are not "
                        "supported together.\n", COLVARS_INPUT_ERROR);
    }
  }

  return COLVARS_OK;
}


int colvarbias_meta::init_well_tempered_params(std::string const &conf)
{
  // for well-tempered metadynamics
  get_keyval(conf, "wellTempered", well_tempered, false);
  get_keyval(conf, "biasTemperature", bias_temperature, -1.0);
  if ((bias_temperature == -1.0) && well_tempered) {
    cvm::error("Error: biasTemperature must be set to a positive value.\n",
               COLVARS_INPUT_ERROR);
  }
  if (well_tempered) {
    cvm::log("Well-tempered metadynamics is used.\n");
    cvm::log("The bias temperature is "+cvm::to_str(bias_temperature)+".\n");
  }
  return COLVARS_OK;
}


int colvarbias_meta::init_ebmeta_params(std::string const &conf)
{
  int error_code = COLVARS_OK;
  // for ebmeta
  target_dist = NULL;
  get_keyval(conf, "ebMeta", ebmeta, false);
  if(ebmeta){
    cvm::main()->cite_feature("Ensemble-biased metadynamics (ebMetaD)");
    if (use_grids && expand_grids) {
      error_code |= cvm::error("Error: expandBoundaries is not supported with "
                               "ebMeta; please allocate wide enough boundaries "
                               "for each colvar ahead of time and set "
                               "targetDistFile accordingly.\n",
                               COLVARS_INPUT_ERROR);
    }
    target_dist = new colvar_grid_scalar();
    error_code |= target_dist->init_from_colvars(colvars);
    std::string target_dist_file;
    get_keyval(conf, "targetDistFile", target_dist_file);
    error_code |= target_dist->read_multicol(target_dist_file,
                                             "ebMeta target histogram");
    cvm::real min_val = target_dist->minimum_value();
    cvm::real max_val = target_dist->maximum_value();
    if (min_val < 0.0) {
      error_code |= cvm::error("Error: Target distribution of EBMetaD "
                               "has negative values!.\n",
                               COLVARS_INPUT_ERROR);
    }
    cvm::real target_dist_min_val;
    get_keyval(conf, "targetDistMinVal", target_dist_min_val, 1/1000000.0);
    if(target_dist_min_val>0 && target_dist_min_val<1){
      target_dist_min_val=max_val*target_dist_min_val;
      target_dist->remove_small_values(target_dist_min_val);
    } else {
      if (target_dist_min_val==0) {
        cvm::log("NOTE: targetDistMinVal is set to zero, the minimum value of the target \n");
        cvm::log(" distribution will be set as the minimum positive value.\n");
        cvm::real min_pos_val = target_dist->minimum_pos_value();
        if (min_pos_val <= 0.0){
          error_code |= cvm::error("Error: Target distribution of EBMetaD has "
                                   "negative or zero minimum positive value.\n",
                                   COLVARS_INPUT_ERROR);
        }
        if (min_val == 0.0){
          cvm::log("WARNING: Target distribution has zero values.\n");
          cvm::log("Zeros will be converted to the minimum positive value.\n");
          target_dist->remove_small_values(min_pos_val);
        }
      } else {
        error_code |= cvm::error("Error: targetDistMinVal must be a value "
                                 "between 0 and 1.\n", COLVARS_INPUT_ERROR);
      }
    }
    // normalize target distribution and multiply by effective volume = exp(differential entropy)
    target_dist->multiply_constant(1.0/target_dist->integral());
    cvm::real volume = cvm::exp(target_dist->entropy());
    target_dist->multiply_constant(volume);
    get_keyval(conf, "ebMetaEquilSteps", ebmeta_equil_steps, ebmeta_equil_steps);
  }

  return error_code;
}

int colvarbias_meta::init_reflection_params(std::string const &conf)
{
  nrefvarsl=0;
  nrefvarsu=0;
  size_t nonpvars;
  size_t i;
  size_t ii=0;
  size_t j;
  size_t jj;

  for ( i = 0; i < num_variables(); i++ ) {
     if (!variables(i)->is_enabled(f_cv_periodic)) {
       ii++;
     }
  }
  nonpvars=ii;

  get_keyval(conf, "useHillsReflection", use_reflection, false);
  if (use_reflection) {

    get_keyval(conf, "reflectionRange", reflection_int, 6.0);
    cvm::log("Reflection range is "+cvm::to_str(reflection_int)+".\n");

    if (get_keyval(conf, "reflectionLowLimitUseCVs", reflection_llimit_cv, reflection_llimit_cv)) {
      nrefvarsl=reflection_llimit_cv.size();
      if(nrefvarsl>num_variables()) cvm::error("Error: number CVs with active lower reflection limit is > num_variables  \n", COLVARS_INPUT_ERROR);
      cvm::log("Using lower limits reflection on "+cvm::to_str(nrefvarsl)+" variables.\n");
    } else {
      nrefvarsl=nonpvars;
      reflection_llimit_cv.resize(nrefvarsl);
      ii=0;
      for (i = 0; i < num_variables(); i++) {
         if (!variables(i)->is_enabled(f_cv_periodic)) {
           reflection_llimit_cv[ii]=i;
           ii++;
         }
      }
      cvm::log("Using all non-periodic variables for lower limits of reflection \n");
    }

    if (reflection_llimit.size()==0) {
      reflection_llimit.resize(nrefvarsl);
    }

    if (get_keyval(conf, "reflectionUpLimitUseCVs", reflection_ulimit_cv, reflection_ulimit_cv)) {
      nrefvarsu=reflection_ulimit_cv.size();
      if(nrefvarsu>num_variables()) cvm::error("Error: number CVs with active upper reflection limit is > num_variables  \n", COLVARS_INPUT_ERROR);
      cvm::log("Using upper limits reflection on "+cvm::to_str(nrefvarsu)+" variables.\n");
    } else {
      nrefvarsu=nonpvars;
      reflection_ulimit_cv.resize(nrefvarsu);
      ii=0;
      for (i = 0; i < num_variables(); i++) {
         if (!variables(i)->is_enabled(f_cv_periodic)) {
           reflection_ulimit_cv[ii]=i;
           ii++;
         }
      }
      cvm::log("Using all non-periodic variables for upper limits of reflection \n");
    }

    if (reflection_ulimit.size()==0) {
      reflection_ulimit.resize(nrefvarsu);
    }

    // use reflection only with scalar variables

    for (i = 0; i < nrefvarsl; i++) {
       if (reflection_llimit_cv[i]>=static_cast<int>(num_variables()) || reflection_llimit_cv[i]<0) {
         cvm::error("Error: CV number is negative or >= num_variables  \n", COLVARS_INPUT_ERROR);
       }
       j=reflection_llimit_cv[i];
       if (variables(j)->value().type()!=colvarvalue::type_scalar) {
         cvm::error("Error: Hills reflection can be used only with scalar variables.\n", COLVARS_INPUT_ERROR);
       }
       if (variables(j)->is_enabled(f_cv_periodic)) {
         cvm::log("Warning: you are using hills reflection with a periodic variable, make sure you are using it far from periodic boundaries \n");
       }
    }

    for (i = 0; i < nrefvarsu; i++) {
       if (reflection_ulimit_cv[i]>=static_cast<int>(num_variables()) || reflection_ulimit_cv[i]<0) {
         cvm::error("Error: CV number is negative or >= num_variables  \n", COLVARS_INPUT_ERROR);
       }
       j=reflection_ulimit_cv[i];
       if (variables(j)->value().type()!=colvarvalue::type_scalar) {
         cvm::error("Error: Hills reflection can be used only with scalar variables.\n", COLVARS_INPUT_ERROR);
       }
       if (variables(j)->is_enabled(f_cv_periodic)) {
         cvm::log("Warning: you are using hills reflection with a periodic variable, make sure you are using it far from periodic boundaries \n");
       }
    }

    // if grids are defined, set by default reflection boundaries as grid boundaries
    if (use_grids) {
      for (i = 0; i < nrefvarsl; i++) {
         ii=reflection_llimit_cv[i];
         reflection_llimit[i]=hills_energy->lower_boundaries[ii].real_value;
      }
      for (i = 0; i < nrefvarsu; i++) {
         ii=reflection_ulimit_cv[i];
         reflection_ulimit[i]=hills_energy->upper_boundaries[ii].real_value;
      }
    }

    if (nrefvarsl>0) {
      if (get_keyval(conf, "reflectionLowLimit", reflection_llimit, reflection_llimit)) {
        for (i = 0; i < nrefvarsl; i++) {
           if (use_grids) {
             ii=reflection_llimit_cv[i];
             cvm:: real bound=hills_energy->lower_boundaries[ii].real_value;
             if (reflection_llimit[i] != bound && reflection_llimit[i] < bound+variables(ii)->width) {
               cvm::error("Error: please set lower reflection limit for CV "+cvm::to_str(ii)+" either at grid lower boundary ("+cvm::to_str(bound)+") or well above it (above "+cvm::to_str(bound+variables(ii)->width)+").\n", COLVARS_INPUT_ERROR);
             }
           }
        }
      } else {
        if (!use_grids) {
          cvm::error("Error: Lower limits for reflection not provided.\n", COLVARS_INPUT_ERROR);
        }
      }
    }

    for (i = 0; i < nrefvarsl; i++) {
       cvm::log("Reflection condition is applied on a lower limit for CV "+cvm::to_str(reflection_llimit_cv[i])+".\n");
       cvm::log("Reflection condition lower limit for this CV is "+cvm::to_str(reflection_llimit[i])+".\n");
    }

    if (nrefvarsu>0) {
      if (get_keyval(conf, "reflectionUpLimit", reflection_ulimit, reflection_ulimit)) {
        for (i = 0; i < nrefvarsu; i++) {
           if (use_grids) {
             ii=reflection_ulimit_cv[i];
             cvm:: real bound=hills_energy->upper_boundaries[ii].real_value;
             if (reflection_ulimit[i] != bound && reflection_ulimit[i] > bound-variables(ii)->width) {
               cvm::error("Error: please set upper reflection limit for CV "+cvm::to_str(ii)+" either at grid upper boundary ("+cvm::to_str(bound)+") or well below it (below "+cvm::to_str(bound-variables(ii)->width)+").\n", COLVARS_INPUT_ERROR);
             }
           }
        }
      } else {
        if (!use_grids) {
          cvm::error("Error: Upper limits for reflection not provided.\n", COLVARS_INPUT_ERROR);
        }
      }
    }

    for (i = 0; i < nrefvarsu; i++) {
       cvm::log("Reflection condition is applied on an upper limit for CV "+cvm::to_str(reflection_ulimit_cv[i])+".\n");
       cvm::log("Reflection condition upper limit for this CV is "+cvm::to_str(reflection_ulimit[i])+".\n");
    }

    // multimensional reflection

    // generate reflection states
    size_t sum;
    sum=1;
    size_t nstates;
    size_t count;

    if (reflection_usel.size()==0) {
      reflection_usel.resize(num_variables(),std::vector<bool>(2));
    }

    if (reflection_l.size()==0) {
      reflection_l.resize(num_variables(),std::vector<cvm::real>(2));
    }

    for (j = 1; j < num_variables(); j++) {
       reflection_usel[j][0]=false;
       reflection_l[j][0]=0.0;
       reflection_usel[j][1]=false;
       reflection_l[j][1]=0.0;
    }

    for (i = 0; i < nrefvarsl; i++) {
       j=reflection_llimit_cv[i];
       reflection_usel[j][0]=true;
       reflection_l[j][0]=reflection_llimit[i];
    }

    for (i = 0; i < nrefvarsu; i++) {
       j=reflection_ulimit_cv[i];
       reflection_usel[j][1]=true;
       reflection_l[j][1]=reflection_ulimit[i];
    }

//  Generate all possible reflection states (e.g. through faces, edges and vertex).
//  Consider for example a cube, the states are:
//  [0,0,1]
//  [0,1,0] [0,1,1]
//  [1,0,0] [1,0,1] [1,1,0] [1,1,1]
//  where 1 means reflect on that coordinate and 0 do not reflect.
//  These states can be generated as:
//  ref_state[0][0]=1
//  ref_state[1][0]=10  ref_state[1][1]=11
//  ref_state[2][0]=100 ref_state[2][1]=101 ref_state[2][2]=110 ref_state[2][3]=111
//  going down along the rows the size ref_state[j].size() is the number of previous states
//  (j-1) plus one.
//  A specific state instead can be generated starting from a power of 10 and then summing
//  the states of the previous rows:
//  ref_state[1][1]=ref_state[1][0]+ref_state[0][0]
//  ref_state[2][1]=ref_state[2][0]+ref_state[0][0]
//  ref_state[2][2]=ref_state[2][0]+ref_state[1][0]
//  ref_state[2][3]=ref_state[2][0]+ref_state[1][1]

    if (ref_state.size()==0) {
      ref_state.resize(num_variables(),std::vector<size_t>(1));
    }
    ref_state[0][0]=1;
    for (j = 1; j < num_variables(); j++) {
      sum*=10;
      nstates=0;
      for (jj = 0; jj < j; jj++) {
            nstates+=ref_state[jj].size();
      }
      nstates++;
      ref_state[j].resize(nstates);
      ref_state[j][0]=sum;
      count=0;
      for (jj = 0; jj < j; jj++) {
         for ( ii = 0; ii < ref_state[jj].size(); ii++) {
            count++;
            ref_state[j][count]=ref_state[j][0]+ref_state[jj][ii];
         }
      }
    }
  }
  return COLVARS_OK;
}

int colvarbias_meta::init_interval_params(std::string const &conf)
{
  use_interval=false;
  nintvarsl=0;
  nintvarsu=0;
  std::vector<int> interval_llimit_cv;
  std::vector<int> interval_ulimit_cv;
  size_t i;
  size_t j;
  // allocate accessory variable of colvar values to be used in energy and forces calculation
  curr_values.resize(num_variables());
  for (i = 0; i < num_variables(); i++) {
    curr_values[i].type(variables(i)->value());
  }

  if (get_keyval(conf, "useHillsInterval", use_interval, use_interval)) {
    if (use_interval) {
      if (get_keyval(conf, "intervalLowLimitUseCVs", interval_llimit_cv, interval_llimit_cv)) {
        nintvarsl=interval_llimit_cv.size();
        cvm::log("Using lower limits interval on "+cvm::to_str(nintvarsl)+" variables.\n");
        for (i = 0; i < nintvarsl; i++) {
           j=interval_llimit_cv[i];
           if (variables(j)->is_enabled(f_cv_periodic)) {
             cvm::log("Warning: you are using interval with a periodic variable, make sure you are using it far from periodic boundaries \n");
           }
        }
      } else {
        cvm::error("Error: Upper limit variables for interval not provided.\n", COLVARS_INPUT_ERROR);
      }
      if (get_keyval(conf, "intervalUpLimitUseCVs", interval_ulimit_cv, interval_ulimit_cv)) {
        nintvarsu=interval_ulimit_cv.size();
        cvm::log("Using upper limits interval on "+cvm::to_str(nintvarsu)+" variables.\n");
        for (i = 0; i < nintvarsu; i++) {
           j=interval_ulimit_cv[i];
           if (variables(j)->is_enabled(f_cv_periodic)) {
             cvm::log("Warning: you are using interval with a periodic variable, make sure you are using it far from periodic boundaries \n");
           }
        }
      } else {
        cvm::error("Error: Upper limits variables for interval not provided.\n", COLVARS_INPUT_ERROR);
      }
      if(nintvarsl>0) {
        if (get_keyval(conf, "intervalLowLimit", interval_llimit, interval_llimit)) {
          for ( i = 0; i < nintvarsl; i++) {
             cvm::log("Hills forces will be removed beyond a lower limit for CV "+cvm::to_str(interval_llimit_cv[i])+".\n");
             cvm::log("Interval condition lower limit for this CV is "+cvm::to_str(interval_llimit[i])+".\n");
          }
        } else {
          cvm::error("Error: Lower limits for interval not provided.\n", COLVARS_INPUT_ERROR);
        }
      }

      if(nintvarsu>0) {
        if (get_keyval(conf, "intervalUpLimit", interval_ulimit, interval_ulimit)) {
          for ( i = 0; i < nintvarsu; i++) {
             cvm::log("Hills forces will be removed beyond an upper limit for CV "+cvm::to_str(interval_ulimit_cv[i])+".\n");
             cvm::log("Interval condition upper limit for this CV is "+cvm::to_str(interval_ulimit[i])+".\n");
          }
        } else {
          cvm::error("Error: Upper limits for interval not provided.\n", COLVARS_INPUT_ERROR);
        }
      }
    }
  } else {
    if (nrefvarsl>0 || nrefvarsu>0) {
      use_interval=true;
      cvm::log("Reflection active: Using by default reflection variables and limits for interval \n");
      nintvarsl=nrefvarsl;
      nintvarsu=nrefvarsu;
      interval_llimit_cv.resize(nintvarsl);
      if (interval_llimit.size()==0) {
        interval_llimit.resize(nintvarsl);
      }
      for (i = 0; i < nintvarsl; i++) {
         interval_llimit_cv[i]=reflection_llimit_cv[i];
         interval_llimit[i]=reflection_llimit[i];
      }
      interval_ulimit_cv.resize(nintvarsu);
      if (interval_ulimit.size()==0) {
        interval_ulimit.resize(nintvarsu);
      }
      for ( i = 0; i < nintvarsu; i++) {
         interval_ulimit_cv[i]=reflection_ulimit_cv[i];
         interval_ulimit[i]=reflection_ulimit[i];
      }
    }
  }

  if (which_int_llimit_cv.size()==0) {
    which_int_llimit_cv.resize(num_variables());
  }
  for ( i = 0; i < num_variables(); i++) {
     which_int_llimit_cv[i]=-1;
  }
  for ( i = 0; i < nintvarsl; i++) {
     j=interval_llimit_cv[i];
     which_int_llimit_cv[j]=i;
  }

  if (which_int_ulimit_cv.size()==0) {
    which_int_ulimit_cv.resize(num_variables());
  }
  for ( i = 0; i < num_variables(); i++) {
     which_int_ulimit_cv[i]=-1;
  }
  for ( i = 0; i < nintvarsu; i++) {
     j=interval_ulimit_cv[i];
     which_int_ulimit_cv[j]=i;
  }
  // use interval only with scalar variables

  for ( i = 0; i < nintvarsl; i++) {
     if (interval_llimit_cv[i]>=static_cast<int>(num_variables()) || interval_llimit_cv[i]<0) {
       cvm::error("Error: CV number is negative or >= num_variables  \n", COLVARS_INPUT_ERROR);
       return COLVARS_INPUT_ERROR;
     }
     j=interval_llimit_cv[i];
     if (variables(j)->value().type()!=colvarvalue::type_scalar) {
       cvm::error("Error: Hills interval can be used only with scalar variables.\n", COLVARS_INPUT_ERROR);
       return COLVARS_INPUT_ERROR;
     }
  }

  for ( i = 0; i < nintvarsu; i++) {
     if (interval_ulimit_cv[i]>=static_cast<int>(num_variables()) || interval_ulimit_cv[i]<0) {
       cvm::error("Error: CV number is negative or >= num_variables  \n", COLVARS_INPUT_ERROR);
       return COLVARS_INPUT_ERROR;
     }
     j=interval_ulimit_cv[i];
     if (variables(j)->value().type()!=colvarvalue::type_scalar) {
       cvm::error("Error: Hills interval can be used only with scalar variables.\n", COLVARS_INPUT_ERROR);
       return COLVARS_INPUT_ERROR;
     }
  }

  return COLVARS_OK;
}


colvarbias_meta::~colvarbias_meta()
{
  colvarbias_meta::clear_state_data();
  colvarproxy *proxy = cvm::main()->proxy;

  proxy->close_output_stream(replica_hills_file);

  proxy->close_output_stream(hills_traj_file_name());

  if (target_dist) {
    delete target_dist;
    target_dist = NULL;
  }
}


int colvarbias_meta::clear_state_data()
{
  if (hills_energy) {
    delete hills_energy;
    hills_energy = NULL;
  }

  if (hills_energy_gradients) {
    delete hills_energy_gradients;
    hills_energy_gradients = NULL;
  }

  hills.clear();
  hills_off_grid.clear();

  return COLVARS_OK;
}


// **********************************************************************
// Hill management member functions
// **********************************************************************

std::list<colvarbias_meta::hill>::const_iterator
colvarbias_meta::add_hill(colvarbias_meta::hill const &h)
{
  hill_iter const hills_end = hills.end();
  hills.push_back(h);
  if (new_hills_begin == hills_end) {
    // if new_hills_begin is unset, set it for the first time
    new_hills_begin = hills.end();
    new_hills_begin--;
  }

  if (use_grids) {

    // also add it to the list of hills that are off-grid, which may
    // need to be computed analytically when the colvar returns
    // off-grid
    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries(h.centers, true);
    if (min_dist < (3.0 * cvm::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(h);
    }
  }

  // output to trajectory (if specified)
  if (b_hills_traj) {
    // Save the current hill to a buffer for further traj output
    hills_traj_os_buf << (hills.back()).output_traj();
  }

  has_data = true;
  return hills.end();
}

bool colvarbias_meta::check_reflection_limits(bool &ah)
{
  size_t i;
  size_t ii;
  for ( i = 0; i < nrefvarsl; i++) {
     ii=reflection_llimit_cv[i];
     if (colvar_values[ii]<reflection_llimit[i]) {
       ah=false;
     }
  }
  for ( i = 0; i < nrefvarsu; i++) {
     ii=reflection_ulimit_cv[i];
     if (colvar_values[ii]>reflection_ulimit[i]) {
       ah=false;
     }
  }
  return ah;
}

int colvarbias_meta::reflect_hill_multid(cvm::real const &h_scale)
{
  size_t i = 0;
  size_t j;
  size_t jj;
  size_t startsum;
  int getsum;
  int check_val;
  size_t numberref;
  size_t startsumk;
  int upordown;
  size_t nkstates;
  size_t kstate;
  size_t k;
  size_t kk;
  size_t countstate;
  bool hill_add;
  int getsumk;
  int checkk;
  size_t state;
  int valk;
  cvm:: real tmps;
  colvarvalue tmp;
  colvarvalue unitary;
  cvm:: real reflection_limit;
  cvm:: real tmpd;

  std::vector<colvarvalue> curr_cv_values(num_variables());
  for (i = 0; i < num_variables(); i++) {
    curr_cv_values[i].type(variables(i)->value());
  }
  for (i = 0; i < num_variables(); i++) {
      curr_cv_values[i] = colvar_values[i];
  }

  // sum over all possible reflection states previously generated,
  // see init

  for ( j = 0; j < num_variables(); j++) {
     startsum=1;
     for (i = 0; i < j; i++) {
        startsum*=10;
     }
     for (jj = 0; jj < ref_state[j].size(); jj++) {
           getsum=startsum;
           check_val=ref_state[j][jj];
           numberref=0;
           startsumk=1;
           for (i = 0; i <= j; i++) {
              upordown=std::floor(check_val/getsum);
              if(check_val-getsum>=0) check_val=check_val-getsum;
              getsum=getsum/10;
              if (upordown==1) {
                numberref++;
                if(numberref>1) startsumk*=10;
              }
           }

           // sum over all possible lower and upper boudary combinations
           // exploiting kstate=ref_state[k][kk]:
           // for just one reflection these are 0(lower boundary) and 1(upper boundary)
           // for two reflections these are 0 1 10 11 (0,0 0,1 1,0 1,1)
           // where 0 is reflect on the two lower boudaries of the two coordinates etc.

           nkstates=2;
           kstate=0;
           for ( k = 0; k <numberref; k++ ) {
              if (k>0)  nkstates=ref_state[k].size();
              for ( kk = 0; kk < nkstates; kk++) {
                 if (k==0 && kk==1) {
                   kstate=1;
                 } else if (k>0) {
                   kstate=ref_state[k][kk];
                 }

                 getsum=startsum;
                 countstate=0;
                 check_val=ref_state[j][jj];
                 hill_add=true;
                 getsumk=startsumk;
                 checkk=kstate;
                 for (i = 0; i <= j; i++) {
                    upordown=std::floor(check_val/getsum);
                    state=num_variables()-1-j+countstate;
                    countstate++;
                    if(check_val-getsum>=0) check_val=check_val-getsum;
                    getsum=getsum/10;
                    if (upordown==1) {
                      tmps=colvar_sigmas[state];
                      tmp=curr_cv_values[state]; // store original current cv value
                      unitary=curr_cv_values[state];
                      unitary.set_ones();
                      valk=std::floor(checkk/getsumk);
                      if(checkk-getsumk>=0) checkk=checkk-getsumk;
                      getsumk=getsumk/10;
                      reflection_limit=reflection_l[state][valk];
                      tmpd=reflection_limit-cvm::real(curr_cv_values[state]);
                      tmpd=std::sqrt(tmpd*tmpd);
                      if (tmpd<reflection_int*tmps && reflection_usel[state][valk] ) { // do mirror within selected range in case upordown=1
                        curr_cv_values[state]=2.0*reflection_limit*unitary-tmp; // reflected cv value
                      } else {
                        hill_add=false;
                      }
                    }
                 }
                 if (hill_add) {
                   switch (comm) {

                   case single_replica:

                     add_hill(hill(cvm::step_absolute(), hill_weight*h_scale, curr_cv_values, colvar_sigmas));

                     break;

                   case multiple_replicas:
                     add_hill(hill(cvm::step_absolute(), hill_weight*h_scale, curr_cv_values, colvar_sigmas, replica_id));
                     std::ostream &replica_hills_os =
                       cvm::proxy->output_stream(replica_hills_file, "replica hills file");
                     if (replica_hills_os) {
                       write_hill(replica_hills_os, hills.back());
                     } else {
                       return cvm::error("Error: in metadynamics bias \""+this->name+"\""+
                                         ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                                         " while writing hills for the other replicas.\n", COLVARS_FILE_ERROR);
                     }
                     break;
                   }

                   for (i = 0; i < num_variables(); i++) {
                      curr_cv_values[i] = colvar_values[i]; // go back to previous values
                   }
                 } else {
                   for (i = 0; i < num_variables(); i++) {
                      curr_cv_values[i] = colvar_values[i]; // go back to previous values
                   }
                 }
              }
           }
     }
  }
  return COLVARS_OK;
}

std::list<colvarbias_meta::hill>::const_iterator
colvarbias_meta::delete_hill(hill_iter &h)
{
  if (cvm::debug()) {
    cvm::log("Deleting hill from the metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
             ", with step number "+
             cvm::to_str(h->it)+(h->replica.size() ?
                                 ", replica id \""+h->replica :
                                 "")+".\n");
  }

  if (use_grids && !hills_off_grid.empty()) {
    for (hill_iter hoff = hills_off_grid.begin();
         hoff != hills_off_grid.end(); hoff++) {
      if (*h == *hoff) {
        hills_off_grid.erase(hoff);
        break;
      }
    }
  }

  if (b_hills_traj) {
    // Save the current hill to a buffer for further traj output
    hills_traj_os_buf << "# DELETED this hill: "
                      << (hills.back()).output_traj()
                      << "\n";
  }

  return hills.erase(h);
}


int colvarbias_meta::update()
{
  int error_code = COLVARS_OK;

  // update base class
  error_code |= colvarbias::update();

  // update the TI estimator (if defined)
  error_code |= colvarbias_ti::update();

  // update grid definition, if needed
  error_code |= update_grid_params();
  // add new biasing energy/forces
  error_code |= update_bias();
  // update grid content to reflect new bias
  error_code |= update_grid_data();

  if (comm != single_replica &&
      (cvm::step_absolute() % replica_update_freq) == 0) {
    // sync with the other replicas (if needed)
    error_code |= replica_share();
  }

  error_code |= calc_energy(NULL);
  error_code |= calc_forces(NULL);

  return error_code;
}


int colvarbias_meta::update_grid_params()
{
  if (use_grids) {

    std::vector<int> curr_bin = hills_energy->get_colvars_index();
    if (cvm::debug()) {
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
               ": current coordinates on the grid: "+
               cvm::to_str(curr_bin)+".\n");
    }

    if (expand_grids) {
      // first of all, expand the grids, if specified
      bool changed_grids = false;
      int const min_buffer =
        (3 * (size_t) cvm::floor(hill_width)) + 1;

      std::vector<int>         new_sizes(hills_energy->sizes());
      std::vector<colvarvalue> new_lower_boundaries(hills_energy->lower_boundaries);
      std::vector<colvarvalue> new_upper_boundaries(hills_energy->upper_boundaries);

      for (size_t i = 0; i < num_variables(); i++) {

        if (! variables(i)->expand_boundaries)
          continue;

        cvm::real &new_lb   = new_lower_boundaries[i].real_value;
        cvm::real &new_ub   = new_upper_boundaries[i].real_value;
        int       &new_size = new_sizes[i];
        bool changed_lb = false, changed_ub = false;

        if (!variables(i)->is_enabled(f_cv_hard_lower_boundary))
          if (curr_bin[i] < min_buffer) {
            int const extra_points = (min_buffer - curr_bin[i]);
            new_lb -= extra_points * variables(i)->width;
            new_size += extra_points;
            // changed offset in this direction => the pointer needs to
            // be changed, too
            curr_bin[i] += extra_points;

            changed_lb = true;
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                     ": new lower boundary for colvar \""+
                     variables(i)->name+"\", at "+
                     cvm::to_str(new_lower_boundaries[i])+".\n");
          }

        if (!variables(i)->is_enabled(f_cv_hard_upper_boundary))
          if (curr_bin[i] > new_size - min_buffer - 1) {
            int const extra_points = (curr_bin[i] - (new_size - 1) + min_buffer);
            new_ub += extra_points * variables(i)->width;
            new_size += extra_points;

            changed_ub = true;
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                     ": new upper boundary for colvar \""+
                     variables(i)->name+"\", at "+
                     cvm::to_str(new_upper_boundaries[i])+".\n");
          }

        if (changed_lb || changed_ub)
          changed_grids = true;
      }

      if (changed_grids) {

        // map everything into new grids

        colvar_grid_scalar *new_hills_energy =
          new colvar_grid_scalar(*hills_energy);
        colvar_grid_gradient *new_hills_energy_gradients =
          new colvar_grid_gradient(*hills_energy_gradients);

        // supply new boundaries to the new grids

        new_hills_energy->lower_boundaries = new_lower_boundaries;
        new_hills_energy->upper_boundaries = new_upper_boundaries;
        new_hills_energy->setup(new_sizes, 0.0, 1);

        new_hills_energy_gradients->lower_boundaries = new_lower_boundaries;
        new_hills_energy_gradients->upper_boundaries = new_upper_boundaries;
        new_hills_energy_gradients->setup(new_sizes, 0.0, num_variables());

        new_hills_energy->map_grid(*hills_energy);
        new_hills_energy_gradients->map_grid(*hills_energy_gradients);

        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy = new_hills_energy;
        hills_energy_gradients = new_hills_energy_gradients;

        curr_bin = hills_energy->get_colvars_index();
        if (cvm::debug())
          cvm::log("Coordinates on the new grid: "+
                   cvm::to_str(curr_bin)+".\n");
      }
    }
  }
  return COLVARS_OK;
}


int colvarbias_meta::update_bias()
{
  colvarproxy *proxy = cvm::main()->proxy;
  // add a new hill if the required time interval has passed
  if (((cvm::step_absolute() % new_hill_freq) == 0) &&
      can_accumulate_data() && is_enabled(f_cvb_history_dependent)) {

    if (cvm::debug()) {
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
               ": adding a new hill at step "+cvm::to_str(cvm::step_absolute())+".\n");
    }

    cvm::real hills_scale=1.0;

    if (ebmeta) {
      hills_scale *= 1.0/target_dist->value(target_dist->get_colvars_index());
      if(cvm::step_absolute() <= ebmeta_equil_steps) {
        cvm::real const hills_lambda =
          (cvm::real(ebmeta_equil_steps - cvm::step_absolute())) /
          (cvm::real(ebmeta_equil_steps));
        hills_scale = hills_lambda + (1-hills_lambda)*hills_scale;
      }
    }

    if (well_tempered) {
      cvm::real hills_energy_sum_here = 0.0;
      if (use_grids) {
        std::vector<int> curr_bin = hills_energy->get_colvars_index();
        hills_energy_sum_here = hills_energy->value(curr_bin);
      } else {
        calc_hills(new_hills_begin, hills.end(), hills_energy_sum_here, NULL);
      }
      hills_scale *= cvm::exp(-1.0*hills_energy_sum_here/(bias_temperature*proxy->boltzmann()));
    }

    // Whether add a hill
    bool add_a_hill=true;

    // Do not add hills beyond reflection borders
    // as just reflected hills must be present
    // beyond those boundaries

    // Check reflection borders: if beyond borders do not add hill

    add_a_hill=check_reflection_limits(add_a_hill);

    if (add_a_hill) {

      switch (comm) {

      case single_replica:

        add_hill(hill(cvm::step_absolute(), hill_weight*hills_scale,
                      colvar_values, colvar_sigmas));

        break;

      case multiple_replicas:
        add_hill(hill(cvm::step_absolute(), hill_weight*hills_scale,
                      colvar_values, colvar_sigmas, replica_id));
        std::ostream &replica_hills_os =
          cvm::proxy->output_stream(replica_hills_file, "replica hills file");
        if (replica_hills_os) {
          write_hill(replica_hills_os, hills.back());
        } else {
          return cvm::error("Error: in metadynamics bias \""+this->name+"\""+
                            ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                            " while writing hills for the other replicas.\n", COLVARS_FILE_ERROR);
        }
        break;
      }

      // add reflected hills if required
      if (use_reflection) {
        reflect_hill_multid(hills_scale);
      }

    }

  }

  return COLVARS_OK;
}


int colvarbias_meta::update_grid_data()
{
  if ((cvm::step_absolute() % grids_freq) == 0) {
    // map the most recent gaussians to the grids
    project_hills(new_hills_begin, hills.end(),
                  hills_energy,    hills_energy_gradients,
                  which_int_llimit_cv, which_int_ulimit_cv,
                  interval_llimit, interval_ulimit);
    new_hills_begin = hills.end();

    // TODO: we may want to condense all into one replicas array,
    // including "this" as the first element
    if (comm == multiple_replicas) {
      for (size_t ir = 0; ir < replicas.size(); ir++) {
        replicas[ir]->project_hills(replicas[ir]->new_hills_begin,
                                    replicas[ir]->hills.end(),
                                    replicas[ir]->hills_energy,
                                    replicas[ir]->hills_energy_gradients,
                                    which_int_llimit_cv, which_int_ulimit_cv,
                                    interval_llimit, interval_ulimit);
        replicas[ir]->new_hills_begin = replicas[ir]->hills.end();
      }
    }
  }

  return COLVARS_OK;
}


int colvarbias_meta::calc_energy(std::vector<colvarvalue> const *values)
{
  size_t ir = 0;

  size_t i;
  int ii;
  std::vector<int> curr_bin;

  curr_values = values ? *values : colvar_values;

  if (use_interval) {
    for (i = 0; i < num_variables(); i++) {
       ii=which_int_llimit_cv[i];
       if (ii>-1 && curr_values[i]<interval_llimit[ii] ) {
         curr_values[i]=interval_llimit[ii];
       }
       ii=which_int_ulimit_cv[i];
       if (ii>-1 && curr_values[i]>interval_ulimit[ii] ) {
         curr_values[i]=interval_ulimit[ii];
       }
    }
    curr_bin = hills_energy->get_colvars_index_bound(curr_values);
  } else {
    curr_bin = hills_energy->get_colvars_index(curr_values);
  }

  for (ir = 0; ir < replicas.size(); ir++) {
    replicas[ir]->bias_energy = 0.0;
  }

  if (hills_energy->index_ok(curr_bin)) {
    // index is within the grid: get the energy from there
    for (ir = 0; ir < replicas.size(); ir++) {

      bias_energy += replicas[ir]->hills_energy->value(curr_bin);
      if (cvm::debug()) {
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                 ": current coordinates on the grid: "+
                 cvm::to_str(curr_bin)+".\n");
        cvm::log("Grid energy = "+cvm::to_str(bias_energy)+".\n");
      }
    }
  } else {
    // off the grid: compute analytically only the hills at the grid's edges
    if (!use_interval) {
      for (ir = 0; ir < replicas.size(); ir++) {
        calc_hills(replicas[ir]->hills_off_grid.begin(),
                   replicas[ir]->hills_off_grid.end(),
                   bias_energy,
                   &curr_values);
      }
    }
  }

  // now include the hills that have not been binned yet (starting
  // from new_hills_begin)

  for (ir = 0; ir < replicas.size(); ir++) {
    calc_hills(replicas[ir]->new_hills_begin,
               replicas[ir]->hills.end(),
               bias_energy,
               &curr_values);
    if (cvm::debug()) {
      cvm::log("Hills energy = "+cvm::to_str(bias_energy)+".\n");
    }
  }

  return COLVARS_OK;
}


int colvarbias_meta::calc_forces(std::vector<colvarvalue> const *values)
{
  size_t ir = 0, ic = 0;
  int ii;
  std::vector<int> curr_bin;
  curr_values = values ? *values : colvar_values;
  std::vector<bool> add_force(num_variables());
  for (ir = 0; ir < replicas.size(); ir++) {
    for (ic = 0; ic < num_variables(); ic++) {
      replicas[ir]->colvar_forces[ic].reset();
    }
  }

  for (ic = 0; ic < num_variables(); ic++) {
    add_force[ic]=true;
    ii=which_int_llimit_cv[ic];
    if (ii>-1) {
      if ( curr_values[ic]<interval_llimit[ii] ) {
        curr_values[ic]=interval_llimit[ii];
        add_force[ic]=false;
      }
    }
    ii=which_int_ulimit_cv[ic];
    if (ii>-1) {
      if ( curr_values[ic]>interval_ulimit[ii] ) {
        add_force[ic]=false;
        curr_values[ic]=interval_ulimit[ii];
      }
    }
  }

  if (use_interval) {
    curr_bin = hills_energy->get_colvars_index_bound(curr_values);
  } else {
    curr_bin = hills_energy->get_colvars_index(curr_values);
  }

  if (hills_energy->index_ok(curr_bin)) {
    for (ir = 0; ir < replicas.size(); ir++) {
      cvm::real const *f = &(replicas[ir]->hills_energy_gradients->value(curr_bin));
      for (ic = 0; ic < num_variables(); ic++) {
        // the gradients are stored, not the forces
        if (add_force[ic]) {
          colvar_forces[ic].real_value += -1.0 * f[ic];
        }
      }
    }
  } else {
    // off the grid: compute analytically only the hills at the grid's edges
    if (!use_interval) {
      for (ir = 0; ir < replicas.size(); ir++) {
        for (ic = 0; ic < num_variables(); ic++) {
           calc_hills_force(ic,
                            replicas[ir]->hills_off_grid.begin(),
                            replicas[ir]->hills_off_grid.end(),
                            colvar_forces,
                            &curr_values);
        }
      }
    }
  }

  // now include the hills that have not been binned yet (starting
  // from new_hills_begin)

  if (cvm::debug()) {
    cvm::log("Metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
             ": adding the forces from the other replicas.\n");
  }

  for (ir = 0; ir < replicas.size(); ir++) {
    for (ic = 0; ic < num_variables(); ic++) {
      if (add_force[ic]) {
        calc_hills_force(ic,
                         replicas[ir]->new_hills_begin,
                         replicas[ir]->hills.end(),
                         colvar_forces,
                         &curr_values);
      }
      if (cvm::debug()) {
        cvm::log("Hills forces = "+cvm::to_str(colvar_forces)+".\n");
      }
    }
  }

  return COLVARS_OK;
}



void colvarbias_meta::calc_hills(colvarbias_meta::hill_iter      h_first,
                                 colvarbias_meta::hill_iter      h_last,
                                 cvm::real                      &energy,
                                 std::vector<colvarvalue> const *values)
{
  size_t i = 0;

  for (hill_iter h = h_first; h != h_last; h++) {

    // compute the gaussian exponent
    cvm::real cv_sqdev = 0.0;
    for (i = 0; i < num_variables(); i++) {
      colvarvalue const &x  = values ? (*values)[i] : colvar_values[i];
      colvarvalue const &center = h->centers[i];
      cvm::real const sigma = h->sigmas[i];
      cv_sqdev += (variables(i)->dist2(x, center)) / (sigma*sigma);
    }

    // compute the gaussian
    if (cv_sqdev > 23.0) {
      // set it to zero if the exponent is more negative than log(1.0E-06)
      h->value(0.0);
    } else {
      h->value(cvm::exp(-0.5*cv_sqdev));
    }
    energy += h->energy();
  }
}


void colvarbias_meta::calc_hills_force(size_t const &i,
                                       colvarbias_meta::hill_iter      h_first,
                                       colvarbias_meta::hill_iter      h_last,
                                       std::vector<colvarvalue>       &forces,
                                       std::vector<colvarvalue> const *values)
{
  // Retrieve the value of the colvar
  colvarvalue const x(values ? (*values)[i] : colvar_values[i]);

  // do the type check only once (all colvarvalues in the hills series
  // were already saved with their types matching those in the
  // colvars)

  hill_iter h;
  switch (x.type()) {

  case colvarvalue::type_scalar:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const sigma = h->sigmas[i];
      forces[i].real_value +=
        ( h->weight() * h->value() * (0.5 / (sigma*sigma)) *
          (variables(i)->dist2_lgrad(x, center)).real_value );
    }
    break;

  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const sigma = h->sigmas[i];
      forces[i].rvector_value +=
        ( h->weight() * h->value() * (0.5 / (sigma*sigma)) *
          (variables(i)->dist2_lgrad(x, center)).rvector_value );
    }
    break;

  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const sigma = h->sigmas[i];
      forces[i].quaternion_value +=
        ( h->weight() * h->value() * (0.5 / (sigma*sigma)) *
          (variables(i)->dist2_lgrad(x, center)).quaternion_value );
    }
    break;

  case colvarvalue::type_vector:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const sigma = h->sigmas[i];
      forces[i].vector1d_value +=
        ( h->weight() * h->value() * (0.5 / (sigma*sigma)) *
          (variables(i)->dist2_lgrad(x, center)).vector1d_value );
    }
    break;

  case colvarvalue::type_notset:
  case colvarvalue::type_all:
  default:
    break;
  }
}


// **********************************************************************
// grid management functions
// **********************************************************************

void colvarbias_meta::project_hills(colvarbias_meta::hill_iter  h_first,
                                    colvarbias_meta::hill_iter  h_last,
                                    colvar_grid_scalar         *he,
                                    colvar_grid_gradient       *hg,
                                    std::vector<int> const &w_int_llimit_cv,
                                    std::vector<int> const &w_int_ulimit_cv,
                                    std::vector<cvm::real> const &int_llimit,
                                    std::vector<cvm::real> const &int_ulimit,
                                    bool print_progress)
{
  if (cvm::debug())
    cvm::log("Metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
             ": projecting hills.\n");

  // TODO: improve it by looping over a small subgrid instead of the whole grid

  std::vector<colvarvalue> new_colvar_values(num_variables());

  std::vector<cvm::real> colvar_forces_scalar(num_variables());

  std::vector<int> he_ix = he->new_index();
  std::vector<int> hg_ix = (hg != NULL) ? hg->new_index() : std::vector<int> (0);
  cvm::real hills_energy_here = 0.0;
  std::vector<colvarvalue> hills_forces_here(num_variables(), 0.0);

  size_t count = 0;
  size_t const print_frequency = ((hills.size() >= 1000000) ? 1 : (1000000/(hills.size()+1)));

  if (hg != NULL) {

    // loop over the points of the grid
    size_t i;
    int ii;
    std::vector<bool> add_force(num_variables());
    for ( ;
          (he->index_ok(he_ix)) && (hg->index_ok(hg_ix));
          count++) {
      for (i = 0; i < num_variables(); i++) {
        add_force[i]=true;
        new_colvar_values[i] = he->bin_to_value_scalar(he_ix[i], i);
        ii=w_int_llimit_cv[i];
        if (ii>-1 ){
          if ( new_colvar_values[i]<int_llimit[ii] ) {
            new_colvar_values[i]=int_llimit[ii];
            add_force[i]=false;
          }
        }
        ii=w_int_ulimit_cv[i];
        if (ii>-1){
          if( new_colvar_values[i]>int_ulimit[ii] ) {
            new_colvar_values[i]=int_ulimit[ii];
            add_force[i]=false;
          }
        }
      }

      // loop over the hills and increment the energy grid locally
      hills_energy_here = 0.0;
      calc_hills(h_first, h_last, hills_energy_here, &new_colvar_values);
      he->acc_value(he_ix, hills_energy_here);

      for (i = 0; i < num_variables(); i++) {
        hills_forces_here[i].reset();
        if (add_force[i]){
          calc_hills_force(i, h_first, h_last, hills_forces_here, &new_colvar_values);
        }
        colvar_forces_scalar[i] = hills_forces_here[i].real_value;
      }
      hg->acc_force(hg_ix, &(colvar_forces_scalar.front()));

      he->incr(he_ix);
      hg->incr(hg_ix);

      if ((count % print_frequency) == 0) {
        if (print_progress) {
          cvm::real const progress = cvm::real(count) / cvm::real(hg->number_of_points());
          std::ostringstream os;
          os.setf(std::ios::fixed, std::ios::floatfield);
          os << std::setw(6) << std::setprecision(2)
             << 100.0 * progress
             << "% done.";
          cvm::log(os.str());
        }
      }
    }

  } else {
    cvm::error("No grid object provided in metadynamics::project_hills()\n",
               COLVARS_BUG_ERROR);
  }

  if (print_progress) {
    cvm::log("100.00% done.\n");
  }

  if (! keep_hills) {
    hills.erase(hills.begin(), hills.end());
  }
}


void colvarbias_meta::recount_hills_off_grid(colvarbias_meta::hill_iter  h_first,
                                             colvarbias_meta::hill_iter  h_last,
                                             colvar_grid_scalar         * /* he */)
{
  hills_off_grid.clear();

  for (hill_iter h = h_first; h != h_last; h++) {
    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries(h->centers, true);
    if (min_dist < (3.0 * cvm::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(*h);
    }
  }
}



// **********************************************************************
// multiple replicas functions
// **********************************************************************


int colvarbias_meta::replica_share()
{
  int error_code = COLVARS_OK;
  // sync with the other replicas (if needed)
  if (comm == multiple_replicas) {
    colvarproxy *proxy = cvm::main()->proxy;
    // reread the replicas registry
    error_code |= update_replicas_registry();
    // empty the output buffer
    error_code |= proxy->flush_output_stream(replica_hills_file);
    error_code |= read_replica_files();
  }
  return error_code;
}


size_t colvarbias_meta::replica_share_freq() const
{
  return replica_update_freq;
}


int colvarbias_meta::update_replicas_registry()
{
  int error_code = COLVARS_OK;

  if (cvm::debug())
    cvm::log("Metadynamics bias \""+this->name+"\""+
             ": updating the list of replicas, currently containing "+
             cvm::to_str(replicas.size())+" elements.\n");

  {
    // copy the whole file into a string for convenience
    std::string line("");
    std::ifstream reg_file(replicas_registry_file.c_str());
    if (reg_file.is_open()) {
      replicas_registry.clear();
      while (colvarparse::getline_nocomments(reg_file, line))
        replicas_registry.append(line+"\n");
    } else {
      error_code |= cvm::error("Error: failed to open file \""+
                               replicas_registry_file+"\" for reading.\n",
                               COLVARS_FILE_ERROR);
    }
  }

  // now parse it
  std::istringstream reg_is(replicas_registry);
  if (reg_is.good()) {

    std::string new_replica("");
    std::string new_replica_file("");
    while ((reg_is >> new_replica) && new_replica.size() &&
           (reg_is >> new_replica_file) && new_replica_file.size()) {

      if (new_replica == this->replica_id) {
        // this is the record for this same replica, skip it
        new_replica_file.clear();
        new_replica.clear();
        continue;
      }

      bool already_loaded = false;
      for (size_t ir = 0; ir < replicas.size(); ir++) {
        if (new_replica == (replicas[ir])->replica_id) {
          // this replica was already added
          if (cvm::debug())
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                     ": skipping a replica already loaded, \""+
                     (replicas[ir])->replica_id+"\".\n");
          already_loaded = true;
          break;
        }
      }

      if (!already_loaded) {
        // add this replica to the registry
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": accessing replica \""+new_replica+"\".\n");
        replicas.push_back(new colvarbias_meta("metadynamics"));
        (replicas.back())->replica_id = new_replica;
        (replicas.back())->replica_list_file = new_replica_file;
        (replicas.back())->replica_state_file = "";
        (replicas.back())->replica_state_file_in_sync = false;

        // Note: the following could become a copy constructor?
        (replicas.back())->name = this->name;
        (replicas.back())->colvars = colvars;
        (replicas.back())->use_grids = use_grids;
        (replicas.back())->dump_fes = false;
        (replicas.back())->expand_grids = false;
        (replicas.back())->rebin_grids = false;
        (replicas.back())->keep_hills = false;
        (replicas.back())->colvar_forces = colvar_forces;

        (replicas.back())->comm = multiple_replicas;

        if (use_grids) {
          (replicas.back())->hills_energy           = new colvar_grid_scalar(colvars);
          (replicas.back())->hills_energy_gradients = new colvar_grid_gradient(colvars);
        }
        if (is_enabled(f_cvb_calc_ti_samples)) {
          (replicas.back())->enable(f_cvb_calc_ti_samples);
          (replicas.back())->colvarbias_ti::init_grids();
        }
        (replicas.back())->update_status = 1;
      }
    }
  } else {
    error_code |= cvm::error("Error: cannot read the replicas registry file \""+
                             replicas_registry+"\".\n", COLVARS_FILE_ERROR);
  }

  // now (re)read the list file of each replica
  for (size_t ir = 0; ir < replicas.size(); ir++) {
    if (cvm::debug())
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ": reading the list file for replica \""+
               (replicas[ir])->replica_id+"\".\n");

    std::ifstream list_is((replicas[ir])->replica_list_file.c_str());
    std::string key;
    std::string new_state_file, new_hills_file;
    if (!(list_is >> key) ||
        !(key == std::string("stateFile")) ||
        !(list_is >> new_state_file) ||
        !(list_is >> key) ||
        !(key == std::string("hillsFile")) ||
        !(list_is >> new_hills_file)) {
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ": failed to read the file \""+
               (replicas[ir])->replica_list_file+"\": will try again after "+
               cvm::to_str(replica_update_freq)+" steps.\n");
      (replicas[ir])->update_status++;
    } else {
      if (new_state_file != (replicas[ir])->replica_state_file) {
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": replica \""+(replicas[ir])->replica_id+
                 "\" has supplied a new state file, \""+new_state_file+
                 "\".\n");
        (replicas[ir])->replica_state_file_in_sync = false;
        (replicas[ir])->replica_state_file = new_state_file;
        (replicas[ir])->replica_hills_file = new_hills_file;
      }
    }
  }

  if (cvm::debug())
    cvm::log("Metadynamics bias \""+this->name+"\": the list of replicas contains "+
             cvm::to_str(replicas.size())+" elements.\n");

  return error_code;
}


int colvarbias_meta::read_replica_files()
{
  // Note: we start from the 2nd replica.
  for (size_t ir = 1; ir < replicas.size(); ir++) {

    // (re)read the state file if necessary
    if ( (! (replicas[ir])->has_data) ||
         (! (replicas[ir])->replica_state_file_in_sync) ) {
      if ((replicas[ir])->replica_state_file.size()) {
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": reading the state of replica \""+
                 (replicas[ir])->replica_id+"\" from file \""+
                 (replicas[ir])->replica_state_file+"\".\n");
        std::ifstream is((replicas[ir])->replica_state_file.c_str());
        if ((replicas[ir])->read_state(is)) {
          // state file has been read successfully
          (replicas[ir])->replica_state_file_in_sync = true;
          (replicas[ir])->update_status = 0;
        } else {
          cvm::log("Failed to read the file \""+
                   (replicas[ir])->replica_state_file+
                   "\": will try again in "+
                   cvm::to_str(replica_update_freq)+" steps.\n");
          (replicas[ir])->replica_state_file_in_sync = false;
          (replicas[ir])->update_status++;
        }
        is.close();
      } else {
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": the state file of replica \""+
                 (replicas[ir])->replica_id+"\" is currently undefined: "
                 "will try again after "+
                 cvm::to_str(replica_update_freq)+" steps.\n");
        (replicas[ir])->update_status++;
      }
    }

    if (! (replicas[ir])->replica_state_file_in_sync) {
      // if a new state file is being read, the hills file is also new
      (replicas[ir])->replica_hills_file_pos = 0;
    }

    // now read the hills added after writing the state file
    if ((replicas[ir])->replica_hills_file.size()) {

      if (cvm::debug())
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": checking for new hills from replica \""+
                 (replicas[ir])->replica_id+"\" in the file \""+
                 (replicas[ir])->replica_hills_file+"\".\n");

      // read hills from the other replicas' files

      std::ifstream is((replicas[ir])->replica_hills_file.c_str());
      if (is.is_open()) {

        // try to resume the previous position (if not the beginning)
        if ((replicas[ir])->replica_hills_file_pos > 0) {
          is.seekg((replicas[ir])->replica_hills_file_pos, std::ios::beg);
        }

        if (!is.is_open()){
          // if fail (the file may have been overwritten), reset this
          // position
          is.clear();
          is.seekg(0, std::ios::beg);
          // reset the counter
          (replicas[ir])->replica_hills_file_pos = 0;
          // schedule to reread the state file
          (replicas[ir])->replica_state_file_in_sync = false;
          // and record the failure
          (replicas[ir])->update_status++;
          cvm::log("Failed to read the file \""+(replicas[ir])->replica_hills_file+
                   "\" at the previous position: will try again in "+
                   cvm::to_str(replica_update_freq)+" steps.\n");
        } else {

          while ((replicas[ir])->read_hill(is)) {
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ": received a hill from replica \""+
                     (replicas[ir])->replica_id+
                     "\" at step "+
                     cvm::to_str(((replicas[ir])->hills.back()).it)+".\n");
          }
          is.clear();
          // store the position for the next read
          (replicas[ir])->replica_hills_file_pos = is.tellg();
          if (cvm::debug()) {
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ": stopped reading file \""+
                     (replicas[ir])->replica_hills_file+
                     "\" at position "+
                     cvm::to_str((replicas[ir])->replica_hills_file_pos)+".\n");
          }

          // test whether this is the end of the file
          is.seekg(0, std::ios::end);
          if (is.tellg() > (replicas[ir])->replica_hills_file_pos + ((std::streampos) 1)) {
            (replicas[ir])->update_status++;
          } else {
            (replicas[ir])->update_status = 0;
          }
        }

      } else {
        cvm::log("Failed to read the file \""+
                 (replicas[ir])->replica_hills_file+
                 "\": will try again in "+
                 cvm::to_str(replica_update_freq)+" steps.\n");
        (replicas[ir])->update_status++;
      }
      is.close();
    }

    size_t const n_flush = (replica_update_freq/new_hill_freq + 1);
    if ((replicas[ir])->update_status > 3*n_flush) {
      // TODO: suspend the calculation?
      cvm::log("WARNING: metadynamics bias \""+this->name+"\""+
               " could not read information from replica \""+
               (replicas[ir])->replica_id+
               "\" after more than "+
               cvm::to_str((replicas[ir])->update_status * replica_update_freq)+
               " steps.  Ensure that it is still running.\n");
    }
  }
  return COLVARS_OK;
}


int colvarbias_meta::set_state_params(std::string const &state_conf)
{
  int error_code = colvarbias::set_state_params(state_conf);

  if (error_code != COLVARS_OK) {
    return error_code;
  }

  colvarparse::get_keyval(state_conf, "keepHills", restart_keep_hills, false,
                          colvarparse::parse_restart);

  if ((!restart_keep_hills) && (cvm::main()->restart_version_number() < 20210604)) {
    if (keep_hills) {
      cvm::log("Warning: could not ensure that keepHills was enabled when "
               "this state file was written; because it is enabled now, "
               "it is assumed that it was also then, but please verify.\n");
      restart_keep_hills = true;
    }
  } else {
    if (restart_keep_hills) {
      cvm::log("This state file/stream contains explicit hills.\n");
    }
  }

  std::string check_replica = "";
  if (colvarparse::get_keyval(state_conf, "replicaID", check_replica,
                              std::string(""), colvarparse::parse_restart) &&
      (check_replica != this->replica_id)) {
    return cvm::error("Error: in the state file , the "
                      "\"metadynamics\" block has a different replicaID ("+
                      check_replica+" instead of "+replica_id+").\n",
                      COLVARS_INPUT_ERROR);
  }

  return COLVARS_OK;
}


template <typename IST, typename GT>
IST & colvarbias_meta::read_grid_data_template_(IST& is, std::string const &key,
                                                GT *grid, GT *backup_grid)
{
  auto const start_pos = is.tellg();
  std::string key_in;
  if (is >> key_in) {
    if ((key != key_in) || !(grid->read_restart(is))) {
      is.clear();
      is.seekg(start_pos);
      is.setstate(std::ios::failbit);
      if (!rebin_grids) {
        if ((backup_grid == nullptr) || (comm == single_replica)) {
          cvm::error("Error: couldn't read grid data for metadynamics bias \""+
                     this->name+"\""+
                     ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                     "; if useGrids was off when the state file was written, "
                     "try enabling rebinGrids now to regenerate the grids.\n", COLVARS_INPUT_ERROR);
        }
      }
    }
  } else {
    is.clear();
    is.seekg(start_pos);
    is.setstate(std::ios::failbit);
  }
  return is;
}


template <typename IST> IST &colvarbias_meta::read_state_data_template_(IST &is)
{
  if (use_grids) {

    colvar_grid_scalar   *hills_energy_backup = NULL;
    colvar_grid_gradient *hills_energy_gradients_backup = NULL;

    if (has_data) {
      if (cvm::debug())
        cvm::log("Backupping grids for metadynamics bias \""+
                 this->name+"\""+
                 ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+".\n");
      hills_energy_backup           = hills_energy;
      hills_energy_gradients_backup = hills_energy_gradients;
      hills_energy                  = new colvar_grid_scalar(colvars);
      hills_energy_gradients        = new colvar_grid_gradient(colvars);
    }

    read_grid_data_template_<IST, colvar_grid_scalar>(is, "hills_energy", hills_energy,
                                                      hills_energy_backup);

    read_grid_data_template_<IST, colvar_grid_gradient>(
        is, "hills_energy_gradients", hills_energy_gradients, hills_energy_gradients_backup);

    if (is) {
      cvm::log("  successfully read the biasing potential and its gradients from grids.\n");
      if (hills_energy_backup != nullptr) {
        // Now that we have successfully updated the grids, delete the backup copies
        delete hills_energy_backup;
        delete hills_energy_gradients_backup;
      }
    } else {
      return is;
    }
  }

  // Save references to the end of the list of existing hills, so that they can
  // be cleared if hills are read successfully from the stream
  bool const existing_hills = !hills.empty();
  size_t const old_hills_size = hills.size();
  hill_iter old_hills_end = hills.end();
  hill_iter old_hills_off_grid_end = hills_off_grid.end();
  if (cvm::debug()) {
    cvm::log("Before reading hills from the state file, there are "+
             cvm::to_str(hills.size())+" hills in memory.\n");
  }

  // Read any hills following the grid data (if any)
  while (read_hill(is)) {
    if (cvm::debug()) {
      cvm::log("Read a previously saved hill under the "
               "metadynamics bias \"" +
               this->name + "\", created at step " + cvm::to_str((hills.back()).it) +
               "; position in stream is " + cvm::to_str(is.tellg()) + ".\n");
    }
  }

  is.clear();

  new_hills_begin = hills.end();
  cvm::log("  successfully read "+cvm::to_str(hills.size() - old_hills_size)+
           " explicit hills from state.\n");

  if (existing_hills) {
    // Prune any hills that pre-existed those just read
    hills.erase(hills.begin(), old_hills_end);
    hills_off_grid.erase(hills_off_grid.begin(), old_hills_off_grid_end);
    if (cvm::debug()) {
      cvm::log("After pruning the old hills, there are now "+
               cvm::to_str(hills.size())+" hills in memory.\n");
    }
  }

  // If rebinGrids is set, rebin the grids based on the current information
  rebin_grids_after_restart();

  if (use_grids) {
    if (!hills_off_grid.empty()) {
      cvm::log(cvm::to_str(hills_off_grid.size())+" hills are near the "
               "grid boundaries: they will be computed analytically "
               "and saved to the state files.\n");
    }
  }

  colvarbias_ti::read_state_data(is);

  if (cvm::debug())
    cvm::log("colvarbias_meta::read_restart() done\n");

  has_data = true;

  if (comm == multiple_replicas) {
    read_replica_files();
  }

  return is;
}


std::istream & colvarbias_meta::read_state_data(std::istream& is)
{
  return read_state_data_template_<std::istream>(is);
}


cvm::memory_stream &colvarbias_meta::read_state_data(cvm::memory_stream &is)
{
  return read_state_data_template_<cvm::memory_stream>(is);
}


void colvarbias_meta::rebin_grids_after_restart()
{
  if (rebin_grids) {

    // allocate new grids (based on the new boundaries and widths just
    // read from the configuration file), and project onto them the
    // grids just read from the restart file

    colvar_grid_scalar   *new_hills_energy =
      new colvar_grid_scalar(colvars);
    colvar_grid_gradient *new_hills_energy_gradients =
      new colvar_grid_gradient(colvars);

    if (cvm::debug()) {
      std::ostringstream tmp_os;
      tmp_os << "hills_energy parameters:\n";
      tmp_os << hills_energy->get_state_params();
      tmp_os << "new_hills_energy parameters:\n";
      tmp_os << new_hills_energy->get_state_params();
      cvm::log(tmp_os.str());
    }

    if (restart_keep_hills && !hills.empty()) {
      // if there are hills, recompute the new grids from them
      cvm::log("Rebinning the energy and forces grids from "+
               cvm::to_str(hills.size())+" hills (this may take a while)...\n");
      project_hills(hills.begin(), hills.end(),
                    new_hills_energy, new_hills_energy_gradients,
                    which_int_llimit_cv, which_int_ulimit_cv,
                    interval_llimit, interval_ulimit, true);
      cvm::log("rebinning done.\n");

    } else {
      // otherwise, use the grids in the restart file
      cvm::log("Rebinning the energy and forces grids "
               "from the grids in the restart file.\n");
      new_hills_energy->map_grid(*hills_energy);
      new_hills_energy_gradients->map_grid(*hills_energy_gradients);
    }

    delete hills_energy;
    delete hills_energy_gradients;
    hills_energy = new_hills_energy;
    hills_energy_gradients = new_hills_energy_gradients;

    // assuming that some boundaries have expanded, eliminate those
    // off-grid hills that aren't necessary any more
    if (!hills.empty())
      recount_hills_off_grid(hills.begin(), hills.end(), hills_energy);
  }
}


template <typename OST>
OST &colvarbias_meta::write_hill_template_(OST &os, colvarbias_meta::hill const &h)
{
  bool const formatted = !std::is_same<OST, cvm::memory_stream>::value;

  if (formatted) {
    os.setf(std::ios::scientific, std::ios::floatfield);
  }

  write_state_data_key(os, "hill", false);

  if (formatted)
    os << "{\n";

  write_state_data_key(os, "step", false);
  if (formatted)
    os << std::setw(cvm::it_width);
  os << h.it;
  if (formatted)
    os << "\n";

  write_state_data_key(os, "weight", false);
  if (formatted)
    os << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width);
  os << h.W;
  if (formatted)
    os << "\n";

  size_t i;
  write_state_data_key(os, "centers", false);
  for (i = 0; i < (h.centers).size(); i++) {
    if (formatted)
      os << " " << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width);
    os << h.centers[i];
  }
  if (formatted)
    os << "\n";

  // For backward compatibility, write the widths instead of the sigmas
  write_state_data_key(os, "widths", false);
  for (i = 0; i < (h.sigmas).size(); i++) {
    if (formatted)
      os << " " << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width);
    os << 2.0 * h.sigmas[i];
  }
  if (formatted)
    os << "\n";

  if (h.replica.size()) {
    write_state_data_key(os, "replicaID", false);
    os << h.replica;
    if (formatted)
      os << "\n";
  }

  if (formatted)
    os << "}\n";

  return os;
}


std::ostream &colvarbias_meta::write_hill(std::ostream &os, colvarbias_meta::hill const &h)
{
  return write_hill_template_<std::ostream>(os, h);
}


cvm::memory_stream &colvarbias_meta::write_hill(cvm::memory_stream &os,
                                                colvarbias_meta::hill const &h)
{
  return write_hill_template_<cvm::memory_stream>(os, h);
}


template <typename IST> IST &hill_stream_error(IST &is, size_t start_pos, std::string const &key)
{
  is.clear();
  is.seekg(start_pos);
  is.setstate(std::ios::failbit);
  cvm::error("Error: in reading data for keyword \"" + key + "\" from stream.\n",
             COLVARS_INPUT_ERROR);
  return is;
}


template <typename IST> IST &colvarbias_meta::read_hill_template_(IST &is)
{
  if (!is)
    return is; // do nothing if failbit is set

  bool const formatted = !std::is_same<IST, cvm::memory_stream>::value;

  auto const start_pos = is.tellg();

  std::string key;
  if (!(is >> key) || (key != "hill")) {
    is.clear();
    is.seekg(start_pos);
    is.setstate(std::ios::failbit);
    return is;
  }

  if (formatted) {
    std::string brace;
    if (!(is >> brace) || (brace != "{")) {
      return hill_stream_error<IST>(is, start_pos, "hill");
    }
  }

  cvm::step_number h_it = 0L;
  cvm::real h_weight = 0.0;
  std::vector<colvarvalue> h_centers(num_variables());
  for (size_t i = 0; i < num_variables(); i++) {
    h_centers[i].type(variables(i)->value());
  }
  std::vector<cvm::real> h_sigmas(num_variables());
  std::string h_replica;

  if (!read_state_data_key(is, "step") || !(is >> h_it)) {
    return hill_stream_error<IST>(is, start_pos, "step");
  }

  if (read_state_data_key(is, "weight")) {
    if (!(is >> h_weight)) {
      return hill_stream_error<IST>(is, start_pos, "weight");
    }
  }

  if (read_state_data_key(is, "centers")) {
    for (size_t i = 0; i < num_variables(); i++) {
      if (!(is >> h_centers[i])) {
        return hill_stream_error<IST>(is, start_pos, "centers");
      }
    }
  }

  if (read_state_data_key(is, "widths")) {
    for (size_t i = 0; i < num_variables(); i++) {
      if (!(is >> h_sigmas[i])) {
        return hill_stream_error<IST>(is, start_pos, "widths");
      }
      // For backward compatibility, read the widths instead of the sigmas
      h_sigmas[i] /= 2.0;
    }
  }

  if (comm != single_replica) {
    if (read_state_data_key(is, "replicaID")) {
      if (!(is >> h_replica)) {
        return hill_stream_error<IST>(is, start_pos, "replicaID");
      }
      if (h_replica != replica_id) {
        cvm::error("Error: trying to read a hill created by replica \"" + h_replica +
                       "\" for replica \"" + replica_id + "\"; did you swap output files?\n",
                   COLVARS_INPUT_ERROR);
        return hill_stream_error<IST>(is, start_pos, "replicaID");
      }
    }
  }

  if (formatted) {
    std::string brace;
    if (!(is >> brace) || (brace != "}")) {
      return hill_stream_error<IST>(is, start_pos, "hill");
    }
  }

  if ((h_it <= state_file_step) && !restart_keep_hills) {
    if (cvm::debug())
      cvm::log("Skipping a hill older than the state file for metadynamics bias \"" + this->name +
               "\"" + ((comm != single_replica) ? ", replica \"" + replica_id + "\"" : "") + "\n");
    return is;
  }

  hill_iter const hills_end = hills.end();
  hills.push_back(hill(h_it, h_weight, h_centers, h_sigmas, h_replica));
  if (new_hills_begin == hills_end) {
    // if new_hills_begin is unset, set it for the first time
    new_hills_begin = hills.end();
    new_hills_begin--;
  }

  if (use_grids) {
    // add this also to the list of hills that are off-grid, which will
    // be computed analytically
    cvm::real const min_dist =
        hills_energy->bin_distance_from_boundaries((hills.back()).centers, true);
    if (min_dist < (3.0 * cvm::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(hills.back());
    }
  }

  has_data = true;
  return is;
}


std::istream &colvarbias_meta::read_hill(std::istream &is)
{
  return read_hill_template_<std::istream>(is);
}


cvm::memory_stream &colvarbias_meta::read_hill(cvm::memory_stream &is)
{
  return read_hill_template_<cvm::memory_stream>(is);
}


int colvarbias_meta::setup_output()
{
  int error_code = COLVARS_OK;

  output_prefix = cvm::output_prefix();
  if (cvm::main()->num_biases_feature(colvardeps::f_cvb_calc_pmf) > 1) {
    // if this is not the only free energy integrator, append
    // this bias's name, to distinguish it from the output of the other
    // biases producing a .pmf file
    output_prefix += ("."+this->name);
  }

  if (comm == multiple_replicas) {

    // TODO: one may want to specify the path manually for intricated filesystems?
    char *pwd = new char[3001];
    if (GETCWD(pwd, 3000) == nullptr) {
      if (pwd != nullptr) { //
        delete[] pwd;
      }
      return cvm::error("Error: cannot get the path of the current working directory.\n",
                        COLVARS_BUG_ERROR);
    }

    replica_list_file =
      (std::string(pwd)+std::string(PATHSEP)+
       this->name+"."+replica_id+".files.txt");
    // replica_hills_file and replica_state_file are those written
    // by the current replica; within the mirror biases, they are
    // those by another replica
    replica_hills_file =
      (std::string(pwd)+std::string(PATHSEP)+
       cvm::output_prefix()+".colvars."+this->name+"."+replica_id+".hills");
    replica_state_file =
      (std::string(pwd)+std::string(PATHSEP)+
       cvm::output_prefix()+".colvars."+this->name+"."+replica_id+".state");
    delete[] pwd;

    // now register this replica

    // first check that it isn't already there
    bool registered_replica = false;
    std::ifstream reg_is(replicas_registry_file.c_str());
    if (reg_is.is_open()) {  // the file may not be there yet
      std::string existing_replica("");
      std::string existing_replica_file("");
      while ((reg_is >> existing_replica) && existing_replica.size() &&
             (reg_is >> existing_replica_file) && existing_replica_file.size()) {
        if (existing_replica == replica_id) {
          // this replica was already registered
          replica_list_file = existing_replica_file;
          reg_is.close();
          registered_replica = true;
          break;
        }
      }
      reg_is.close();
    }

    // if this replica was not included yet, we should generate a
    // new record for it: but first, we write this replica's files,
    // for the others to read

    // open the "hills" buffer file
    reopen_replica_buffer_file();

    // write the state file (so that there is always one available)
    write_replica_state_file();

    // schedule to read the state files of the other replicas
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      (replicas[ir])->replica_state_file_in_sync = false;
    }

    // if we're running without grids, use a growing list of "hills" files
    // otherwise, just one state file and one "hills" file as buffer
    std::ostream &list_os = cvm::proxy->output_stream(replica_list_file, "replica list file");
    if (list_os) {
      list_os << "stateFile " << replica_state_file << "\n";
      list_os << "hillsFile " << replica_hills_file << "\n";
      cvm::proxy->close_output_stream(replica_list_file);
    } else {
      error_code |= COLVARS_FILE_ERROR;
    }

    // finally, add a new record for this replica to the registry
    if (! registered_replica) {
      std::ofstream reg_os(replicas_registry_file.c_str(), std::ios::app);
      if (!reg_os) {
        return cvm::get_error();
      }
      reg_os << replica_id << " " << replica_list_file << "\n";
      cvm::proxy->close_output_stream(replicas_registry_file);
    }
  }

  if (b_hills_traj) {
    std::ostream &hills_traj_os =
      cvm::proxy->output_stream(hills_traj_file_name(), "hills trajectory file");
    if (!hills_traj_os) {
      error_code |= COLVARS_FILE_ERROR;
    }
  }

  return error_code;
}


std::string const colvarbias_meta::hills_traj_file_name() const
{
  return std::string(cvm::output_prefix()+
                     ".colvars."+this->name+
                     ( (comm != single_replica) ?
                       ("."+replica_id) :
                       ("") )+
                     ".hills.traj");
}


std::string const colvarbias_meta::get_state_params() const
{
  std::ostringstream os;
  if (keep_hills) {
    os << "keepHills on" << "\n";
  }
  if (this->comm != single_replica) {
    os << "replicaID " << this->replica_id << "\n";
  }
  return (colvarbias::get_state_params() + os.str());
}


template <typename OST> OST &colvarbias_meta::write_state_data_template_(OST &os)
{
  if (use_grids) {

    // this is a very good time to project hills, if you haven't done
    // it already!
    project_hills(new_hills_begin, hills.end(),
                  hills_energy,    hills_energy_gradients,
                  which_int_llimit_cv, which_int_ulimit_cv,
                  interval_llimit, interval_ulimit);
    new_hills_begin = hills.end();

    // write down the grids to the restart file
    write_state_data_key(os, "hills_energy");
    hills_energy->write_restart(os);
    write_state_data_key(os, "hills_energy_gradients");
    hills_energy_gradients->write_restart(os);
  }

  if ((!use_grids) || keep_hills) {
    // write all hills currently in memory
    for (std::list<hill>::const_iterator h = this->hills.begin(); h != this->hills.end(); h++) {
      write_hill(os, *h);
    }
  } else {
    // write just those that are near the grid boundaries
    for (std::list<hill>::const_iterator h = this->hills_off_grid.begin();
         h != this->hills_off_grid.end(); h++) {
      write_hill(os, *h);
    }
  }

  colvarbias_ti::write_state_data(os);
  return os;
}


std::ostream & colvarbias_meta::write_state_data(std::ostream& os)
{
  return write_state_data_template_<std::ostream>(os);
}


cvm::memory_stream &colvarbias_meta::write_state_data(cvm::memory_stream &os)
{
  return write_state_data_template_<cvm::memory_stream>(os);
}


int colvarbias_meta::write_state_to_replicas()
{
  int error_code = COLVARS_OK;
  if (comm != single_replica) {
    error_code |= write_replica_state_file();
    error_code |= reopen_replica_buffer_file();
    // schedule to reread the state files of the other replicas
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      (replicas[ir])->replica_state_file_in_sync = false;
    }
  }
  return error_code;
}


int colvarbias_meta::write_output_files()
{
  colvarbias_ti::write_output_files();
  if (dump_fes) {
    write_pmf();
  }
  if (b_hills_traj) {
    std::ostream &hills_traj_os =
        cvm::proxy->output_stream(hills_traj_file_name(), "hills trajectory file");
    hills_traj_os << hills_traj_os_buf.str();
    cvm::proxy->flush_output_stream(hills_traj_file_name());
    // clear the buffer
    hills_traj_os_buf.str("");
    hills_traj_os_buf.clear();
  }
  return COLVARS_OK;
}


void colvarbias_meta::write_pmf()
{
  colvarproxy *proxy = cvm::main()->proxy;
  // allocate a new grid to store the pmf
  colvar_grid_scalar *pmf = new colvar_grid_scalar(*hills_energy);
  pmf->setup();

  if ((comm == single_replica) || (dump_replica_fes)) {
    // output the PMF from this instance or replica
    pmf->reset();
    pmf->add_grid(*hills_energy);

    if (ebmeta) {
      int nt_points=pmf->number_of_points();
      for (int i = 0; i < nt_points; i++) {
         cvm::real pmf_val=0.0;
         cvm::real target_val=target_dist->value(i);
         if (target_val>0) {
           pmf_val=pmf->value(i);
           pmf_val=pmf_val + proxy->target_temperature() * proxy->boltzmann() * cvm::logn(target_val);
         }
         pmf->set_value(i,pmf_val);
      }
    }

    cvm::real const max = pmf->maximum_value();
    pmf->add_constant(-1.0 * max);
    pmf->multiply_constant(-1.0);
    if (well_tempered) {
      cvm::real const well_temper_scale = (bias_temperature + proxy->target_temperature()) / bias_temperature;
      pmf->multiply_constant(well_temper_scale);
    }
    {
      std::string const fes_file_name(this->output_prefix +
                                      ((comm != single_replica) ? ".partial" : "") +
                                      (dump_fes_save ?
                                       "."+cvm::to_str(cvm::step_absolute()) : "") +
                                      ".pmf");
      pmf->write_multicol(fes_file_name, "PMF file");
    }
  }

  if (comm != single_replica) {
    // output the combined PMF from all replicas
    pmf->reset();
    // current replica already included in the pools of replicas
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      pmf->add_grid(*(replicas[ir]->hills_energy));
    }

    if (ebmeta) {
      int nt_points=pmf->number_of_points();
      for (int i = 0; i < nt_points; i++) {
         cvm::real pmf_val=0.0;
         cvm::real target_val=target_dist->value(i);
         if (target_val>0) {
           pmf_val=pmf->value(i);
           pmf_val=pmf_val + proxy->target_temperature() * proxy->boltzmann() * cvm::logn(target_val);
         }
         pmf->set_value(i,pmf_val);
      }
    }

    cvm::real const max = pmf->maximum_value();
    pmf->add_constant(-1.0 * max);
    pmf->multiply_constant(-1.0);
    if (well_tempered) {
      cvm::real const well_temper_scale = (bias_temperature + proxy->target_temperature()) / bias_temperature;
      pmf->multiply_constant(well_temper_scale);
    }
    std::string const fes_file_name(this->output_prefix +
                                    (dump_fes_save ?
                                     "."+cvm::to_str(cvm::step_absolute()) : "") +
                                    ".pmf");
    pmf->write_multicol(fes_file_name, "partial PMF file");
  }

  delete pmf;
}



int colvarbias_meta::write_replica_state_file()
{
  colvarproxy *proxy = cvm::proxy;

  if (cvm::debug()) {
    cvm::log("Writing replica state file for bias \""+name+"\"\n");
  }

  int error_code = COLVARS_OK;

  // Write to temporary state file
  std::string const tmp_state_file(replica_state_file+".tmp");
  error_code |= proxy->remove_file(tmp_state_file);
  std::ostream &rep_state_os = cvm::proxy->output_stream(tmp_state_file, "temporary state file");
  if (rep_state_os) {
    if (!write_state(rep_state_os)) {
      error_code |= cvm::error("Error: in writing to temporary file \""+
                               tmp_state_file+"\".\n", COLVARS_FILE_ERROR);
    }
  }
  error_code |= proxy->close_output_stream(tmp_state_file);

  error_code |= proxy->rename_file(tmp_state_file, replica_state_file);

  return error_code;
}


int colvarbias_meta::reopen_replica_buffer_file()
{
  int error_code = COLVARS_OK;
  colvarproxy *proxy = cvm::proxy;
  if (proxy->output_stream(replica_hills_file, "replica hills file")) {
    error_code |= proxy->close_output_stream(replica_hills_file);
  }
  error_code |= proxy->remove_file(replica_hills_file);
  std::ostream &replica_hills_os = proxy->output_stream(replica_hills_file, "replica hills file");
  if (replica_hills_os) {
    replica_hills_os.setf(std::ios::scientific, std::ios::floatfield);
  } else {
    error_code |= COLVARS_FILE_ERROR;
  }
  return error_code;
}


std::string colvarbias_meta::hill::output_traj()
{
  std::ostringstream os;
  os.setf(std::ios::fixed, std::ios::floatfield);
  os << std::setw(cvm::it_width) << it << " ";

  os.setf(std::ios::scientific, std::ios::floatfield);

  size_t i;
  os << "  ";
  for (i = 0; i < centers.size(); i++) {
    os << " ";
    os << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)  << centers[i];
  }

  os << "  ";
  for (i = 0; i < sigmas.size(); i++) {
    os << " ";
    os << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width) << sigmas[i];
  }

  os << "  ";
  os << std::setprecision(cvm::en_prec)
     << std::setw(cvm::en_width) << W << "\n";

  return os.str();
}


colvarbias_meta::hill::hill(cvm::step_number it_in,
                            cvm::real W_in,
                            std::vector<colvarvalue> const &cv_values,
                            std::vector<cvm::real> const &cv_sigmas,
                            std::string const &replica_in)
  : it(it_in),
    sW(1.0),
    W(W_in),
    centers(cv_values.size()),
    sigmas(cv_values.size()),
    replica(replica_in)
{
  hill_value = 0.0;
  for (size_t i = 0; i < cv_values.size(); i++) {
    centers[i].type(cv_values[i]);
    centers[i] = cv_values[i];
    sigmas[i] = cv_sigmas[i];
  }
  if (cvm::debug()) {
    cvm::log("New hill, applied to "+cvm::to_str(cv_values.size())+
             " collective variables, with centers "+
             cvm::to_str(centers)+", sigmas "+
             cvm::to_str(sigmas)+" and weight "+
             cvm::to_str(W)+".\n");
  }
}


colvarbias_meta::hill::hill(colvarbias_meta::hill const &h)
  : it(h.it),
    hill_value(0.0),
    sW(1.0),
    W(h.W),
    centers(h.centers),
    sigmas(h.sigmas),
    replica(h.replica)
{
  hill_value = 0.0;
}


colvarbias_meta::hill &
colvarbias_meta::hill::operator = (colvarbias_meta::hill const &h)
{
  it = h.it;
  hill_value = 0.0;
  sW = 1.0;
  W = h.W;
  centers = h.centers;
  sigmas = h.sigmas;
  replica = h.replica;
  hill_value = h.hill_value;
  return *this;
}


colvarbias_meta::hill::~hill()
{}
