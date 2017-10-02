// -*- c++ -*-

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

// used to set the absolute path of a replica file
#if defined(WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define CHDIR ::_chdir
#define GETCWD ::_getcwd
#define PATHSEP "\\"
#else
#include <unistd.h>
#define CHDIR ::chdir
#define GETCWD ::getcwd
#define PATHSEP "/"
#endif


#include "colvar.h"
#include "colvarbias_meta.h"


colvarbias_meta::colvarbias_meta(char const *key)
  : colvarbias(key)
{
  new_hills_begin = hills.end();
  hills_traj_os = NULL;
  replica_hills_os = NULL;
}


int colvarbias_meta::init(std::string const &conf)
{
  colvarbias::init(conf);

  enable(f_cvb_calc_pmf);

  get_keyval(conf, "hillWeight", hill_weight, 0.0);
  if (hill_weight > 0.0) {
    enable(f_cvb_apply_force);
  } else {
    cvm::error("Error: hillWeight must be provided, and a positive number.\n", INPUT_ERROR);
  }

  get_keyval(conf, "newHillFrequency", new_hill_freq, 1000);
  if (new_hill_freq > 0) {
    enable(f_cvb_history_dependent);
  }

  get_keyval(conf, "hillWidth", hill_width, std::sqrt(2.0 * PI) / 2.0);
  cvm::log("Half-widths of the Gaussian hills (sigma's):\n");
  for (size_t i = 0; i < num_variables(); i++) {
    cvm::log(variables(i)->name+std::string(": ")+
             cvm::to_str(0.5 * variables(i)->width * hill_width));
  }

  {
    bool b_replicas = false;
    get_keyval(conf, "multipleReplicas", b_replicas, false);
    if (b_replicas)
      comm = multiple_replicas;
    else
      comm = single_replica;
  }

  // in all cases, the first replica is this bias itself
  if (replicas.size() == 0) {
    replicas.push_back(this);
  }

  get_keyval(conf, "useGrids", use_grids, true);

  if (use_grids) {
    get_keyval(conf, "gridsUpdateFrequency", grids_freq, new_hill_freq);
    get_keyval(conf, "rebinGrids", rebin_grids, false);

    expand_grids = false;
    size_t i;
    for (i = 0; i < num_variables(); i++) {
      variables(i)->enable(f_cv_grid);
      if (variables(i)->expand_boundaries) {
        expand_grids = true;
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                 ": Will expand grids when the colvar \""+
                 variables(i)->name+"\" approaches its boundaries.\n");
      }
    }

    get_keyval(conf, "keepHills", keep_hills, false);
    if (! get_keyval(conf, "writeFreeEnergyFile", dump_fes, true))
      get_keyval(conf, "dumpFreeEnergyFile", dump_fes, true, colvarparse::parse_silent);
    if (get_keyval(conf, "saveFreeEnergyFile", dump_fes_save, false, colvarparse::parse_silent)) {
      cvm::log("Option \"saveFreeEnergyFile\" is deprecated, "
               "please use \"keepFreeEnergyFile\" instead.");
    }
    get_keyval(conf, "keepFreeEnergyFiles", dump_fes_save, dump_fes_save);

    hills_energy           = new colvar_grid_scalar(colvars);
    hills_energy_gradients = new colvar_grid_gradient(colvars);
  } else {
    rebin_grids = false;
    keep_hills = false;
    dump_fes = false;
    dump_fes_save = false;
    dump_replica_fes = false;

    hills_energy           = NULL;
    hills_energy_gradients = NULL;
  }

  if (comm != single_replica) {

    if (expand_grids)
      cvm::fatal_error("Error: expandBoundaries is not supported when "
                       "using more than one replicas; please allocate "
                       "wide enough boundaries for each colvar"
                       "ahead of time.\n");

    if (get_keyval(conf, "dumpPartialFreeEnergyFile", dump_replica_fes, false)) {
      if (dump_replica_fes && (! dump_fes)) {
        cvm::log("Enabling \"dumpFreeEnergyFile\".\n");
      }
    }

    get_keyval(conf, "replicaID", replica_id, std::string(""));
    if (!replica_id.size())
      cvm::error("Error: replicaID must be defined "
                 "when using more than one replica.\n", INPUT_ERROR);

    get_keyval(conf, "replicasRegistry",
               replicas_registry_file,
               (this->name+".replicas.registry.txt"));

    get_keyval(conf, "replicaUpdateFrequency",
               replica_update_freq, new_hill_freq);

    if (keep_hills)
      cvm::log("Warning: in metadynamics bias \""+this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
               ": keepHills with more than one replica can lead to a very "
               "large amount of input/output and slow down your calculations.  "
               "Please consider disabling it.\n");

  }

  get_keyval(conf, "writeHillsTrajectory", b_hills_traj, false);

  init_well_tempered_params(conf);

  // read hills inversion or reflection

  init_inversion_params(conf);
  init_reflection_params(conf);
  init_interval_params(conf);

  // init ebmeta and hills scale kernel
  init_ebmeta_params(conf);
  init_kernel_params(conf);

  if (cvm::debug())
    cvm::log("Done initializing the metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+".\n");

  return COLVARS_OK;
}

int colvarbias_meta::init_kernel_params(std::string const &conf)
{
  // use specified kernel to scale hills
  kernel_coupling_time = 0.0;
  kernel_type = kt_none;
  default_kernel_ebmeta=false;
  std::string kernel_type_str;
  get_keyval(conf, "hillsKernel", scale_kernel, false);
  get_keyval(conf, "hillsKernelType", kernel_type_str, to_lower_cppstr(std::string("inverseSqrtTime")));
  kernel_type_str = to_lower_cppstr(kernel_type_str);
  if (kernel_type_str == to_lower_cppstr(std::string("inverseSqrtTime"))) {
    kernel_type = kt_inv_sqrt_time;
  } else if (kernel_type_str == to_lower_cppstr(std::string("none"))) {
    kernel_type = kt_none;
  }

  if (scale_kernel) {
    cvm::log("A scaling time kernel for the hills of metadynamics is used.\n");
    if (well_tempered) {
      cvm::log("WARNING: you are using a scaling kernel for the hills of metadynamics \n");
      cvm::log("together with Well-tempered metadynamics; are you sure this is a good idea?\n");
    }
    switch (kernel_type) {
    case kt_none:
      cvm::error("Error: undefined kernel type.\n", INPUT_ERROR);
      return INPUT_ERROR;
      break;
    case kt_inv_sqrt_time:
    case kt_ntot:
      provide(f_cvb_history_dependent);
      break;
    }
  }

  get_keyval(conf, "hillsKernelCouplingTime", kernel_coupling_time, 0.0);
  if (scale_kernel && kernel_coupling_time <= 0.0) {
    if (ebmeta) {
      default_kernel_ebmeta=true;
      cvm::log("Default ebmeta hillsKernelCouplingTime is used.\n");
    } else {
      cvm::error("Error: a positive couplingTime must be provided.\n", INPUT_ERROR);
      return INPUT_ERROR;
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
    cvm::fatal_error("Error: biasTemperature is not set.\n");
  }
  if (well_tempered) {
    cvm::log("Well-tempered metadynamics is used.\n");
    cvm::log("The bias temperature is "+cvm::to_str(bias_temperature)+".\n");
  }
  return COLVARS_OK;
}


int colvarbias_meta::init_ebmeta_params(std::string const &conf)
{
  // for ebmeta
  // calculate gaussian scale factor in units of hills steps
  gauss_factor = 1.0;
  for (size_t i = 0; i < num_variables(); i++) {
      cvm:: real sigma=0.5*variables(i)->width*hill_width;
      gauss_factor*=sigma*std::sqrt(2.0 * PI);
  } 
  gauss_factor = cvm::kt()/(hill_weight*gauss_factor);
  
  target_dist = NULL;
  target_error = NULL;
  get_keyval(conf, "ebMeta", ebmeta, false);
  if(ebmeta){
    if (use_grids && expand_grids) {
      cvm::fatal_error("Error: expandBoundaries is not supported with "
                       "ebMeta please allocate wide enough boundaries for "
                       "each colvar ahead of time and set targetdistfile "
                       "accordingly. \n");
    }
    target_dist = new colvar_grid_scalar();
    target_dist->init_from_colvars(colvars);
    if (get_keyval(conf, "targetDistFile", target_dist_file)) {
      std::ifstream targetdiststream(target_dist_file.c_str());
      target_dist->read_multicol(targetdiststream);
    } else {
      cvm::log("NOTE: targetDistFile not found; using uniform target distribution by default .\n");
      int nt_points=target_dist->number_of_points();
      for (size_t i = 0; i < nt_points; i++) { 
         target_dist->set_array_value(i,1.0);
      }
    }
    cvm::real min_val = target_dist->minimum_value();
    if(min_val<0){
      cvm::error("Error: Target distribution of ebMeta "
                 "has negative values!.\n", INPUT_ERROR);
    }
    get_keyval(conf, "ebmetaLowLimitNCVs", nebmvarsl, 0);
    get_keyval(conf, "ebmetaUpLimitNCVs", nebmvarsu, 0);

    if (nebmvarsl>0) {
      if (ebmeta_llimit_cv.size()==0) {
        ebmeta_llimit_cv.resize(nebmvarsl);
      }
      if (get_keyval(conf, "ebmetaLowLimitUseCVs", ebmeta_llimit_cv, ebmeta_llimit_cv)) {
        if (ebmeta_llimit.size()==0) {
          ebmeta_llimit.resize(nebmvarsl);
        }
      } else {
        cvm::error("Error in ebmeta input: which CVs have a lower limit not provided.\n", INPUT_ERROR);
        return INPUT_ERROR;        
      }
      if (get_keyval(conf, "ebmetaLowLimit", ebmeta_llimit, ebmeta_llimit)) {
        for (size_t i = 0; i < nebmvarsl; i++) {
           if (ebmeta_llimit_cv[i]>=num_variables() || ebmeta_llimit_cv[i]<0) {
             cvm::error("Error: CV number is negative or >= num_variables  \n", INPUT_ERROR);
             return INPUT_ERROR;
           }
           cvm::log("ebmeta is applied with a lower limit for CV "+cvm::to_str(ebmeta_llimit_cv[i])+".\n");
           cvm::log("ebmeta lower limit for this CV is "+cvm::to_str(ebmeta_llimit[i])+".\n");
        }
      } else {
        cvm::error("Error: Lower limits for ebmeta CVs not provided.\n", INPUT_ERROR);
        return INPUT_ERROR;
      }
    }

    if (nebmvarsu>0) {
      if (ebmeta_ulimit_cv.size()==0) {
        ebmeta_ulimit_cv.resize(nebmvarsu);
      }
      if (get_keyval(conf, "ebmetaUpLimitUseCVs", ebmeta_ulimit_cv, ebmeta_ulimit_cv)) {
        if (ebmeta_ulimit.size()==0) {
          ebmeta_ulimit.resize(nebmvarsu);
        }
      } else {
        cvm::error("Error in ebmeta input: which CVs have a upper limit not provided.\n", INPUT_ERROR);
        return INPUT_ERROR;
      }
      if (get_keyval(conf, "ebmetaUpLimit", ebmeta_ulimit, ebmeta_ulimit)) {
        for (size_t i = 0; i < nebmvarsu; i++) {
           if (ebmeta_ulimit_cv[i]>=num_variables() || ebmeta_ulimit_cv[i]<0) {
             cvm::error("Error: CV number is negative or >= num_variables  \n", INPUT_ERROR);
             return INPUT_ERROR;
           }
           cvm::log("ebmeta is applied with a upper limit for CV "+cvm::to_str(ebmeta_ulimit_cv[i])+".\n");
           cvm::log("ebmeta upper limit for this CV is "+cvm::to_str(ebmeta_ulimit[i])+".\n");
        }
      } else {
        cvm::error("Error: Upper limits for ebmeta CVs not provided.\n", INPUT_ERROR);
        return INPUT_ERROR;
      }
    }

    if (nebmvarsl>0 && nrefvarsl>0 ){
      cvm::log("NOTE: please verify that ebmeta and reflection lower limits are the same for the same CV \n");      
    }

    if (nebmvarsu>0 && nrefvarsu>0 ){
      cvm::log("NOTE: please verify that ebmeta and reflection upper limits are the same for the same CV \n");
    }

    // filter target distribution according to selected limits

    for (size_t j = 0; j < nebmvarsl; j++) {
      size_t i=ebmeta_llimit_cv[j];
      for (std::vector<int> ix = target_dist->new_index(); target_dist->index_ok(ix); target_dist->incr(ix) ) {
         cvm::real cv_value=target_dist->bin_to_value_scalar(ix[i], i);
         if (cv_value < ebmeta_llimit[j]){
           target_dist->set_value(ix, 0.0, 0);
         }
      }
    }

    for (size_t j = 0; j < nebmvarsu; j++) {
      size_t i=ebmeta_ulimit_cv[j];
      for (std::vector<int> ix = target_dist->new_index(); target_dist->index_ok(ix); target_dist->incr(ix) ) {
         cvm::real cv_value=target_dist->bin_to_value_scalar(ix[i], i);
         if (cv_value > ebmeta_ulimit[j]){
           target_dist->set_value(ix, 0.0, 0);
         }
      }
    }

    // filter target distribution according to reflection boundaries

    for (size_t j = 0; j < nrefvarsl; j++) {
      size_t i=reflection_llimit_cv[j];
      for (std::vector<int> ix = target_dist->new_index(); target_dist->index_ok(ix); target_dist->incr(ix) ) {
         cvm::real cv_value=target_dist->bin_to_value_scalar(ix[i], i);
         if (cv_value < reflection_llimit[j]){
           target_dist->set_value(ix, 0.0, 0); 
         } 
      }
    }

    for (size_t j = 0; j < nrefvarsu; j++) {
      size_t i=reflection_ulimit_cv[j];
      for (std::vector<int> ix = target_dist->new_index(); target_dist->index_ok(ix); target_dist->incr(ix) ) {
         cvm::real cv_value=target_dist->bin_to_value_scalar(ix[i], i);
         if (cv_value > reflection_ulimit[j]){
           target_dist->set_value(ix, 0.0, 0);
         }
      }
    }

    // normalize target distribution
    cvm::real intval = 1.0/target_dist->integral();
    target_dist->multiply_constant(intval);
    get_keyval(conf, "ebMetaEquilSteps", ebmeta_equil_steps, 0);
    get_keyval(conf, "ebMetaError", ebmetaerror, false);
    if(ebmetaerror){
      ebmeta_out=false;
      if(get_keyval(conf, "ebmetaOutFile", ebmeta_out_file)){
        ebmetaoutfile.open (ebmeta_out_file);
        ebmeta_out=true;
      } 
      ebmeta_nconst=0.0;
      get_keyval(conf, "ebMetaFixErrorBound", ebmeta_fix_bound, 1.5);
      int const npoints=target_dist->number_of_points();
      target_error = new colvar_grid_scalar();
      target_error->init_from_colvars(colvars);
      if (get_keyval(conf, "targetErrorFile", target_error_file)) {
        std::ifstream targeterrorstream(target_error_file.c_str());
        target_error->read_multicol(targeterrorstream);
      } else {
        cvm::log("NOTE: targetErrorFile not found; using 0.25*target_dist as default error .\n");
        int nt_points=target_error->number_of_points();
        for (size_t i = 0; i < nt_points; i++) {
           cvm:: real error_val=0.25*target_dist->array_value(i)/intval; 
           target_error->set_array_value(i,error_val);
        }
      }
      
      cvm::real error_min = target_error->minimum_value();
      if(error_min<0){
        cvm::error("Error: Target error of ebMeta "
                   "has negative values!.\n", INPUT_ERROR);

      }

      // filter target error according to selected limits
      for (size_t j = 0; j < nebmvarsl; j++) {
        size_t i=ebmeta_llimit_cv[j];
        for (std::vector<int> ix = target_error->new_index(); target_error->index_ok(ix); target_error->incr(ix) ) {
           cvm::real cv_value=target_error->bin_to_value_scalar(ix[i], i);
           if (cv_value < ebmeta_llimit[j]){
             target_error->set_value(ix, 0.0, 0);
           }
        }
      }
   
      for (size_t j = 0; j < nebmvarsu; j++) {
        size_t i=ebmeta_ulimit_cv[j];
        for (std::vector<int> ix = target_error->new_index(); target_error->index_ok(ix); target_error->incr(ix) ) {
           cvm::real cv_value=target_error->bin_to_value_scalar(ix[i], i);
           if (cv_value > ebmeta_ulimit[j]){
             target_error->set_value(ix , 0.0, 0);
           }
        }
      }

      // filter target error according to reflection boundaries
      for (size_t j = 0; j < nrefvarsl; j++) {
        size_t i=reflection_llimit_cv[j];
        for (std::vector<int> ix = target_error->new_index(); target_error->index_ok(ix); target_error->incr(ix) ) {
           cvm::real cv_value=target_error->bin_to_value_scalar(ix[i], i);
           if (cv_value < reflection_llimit[j]){
             target_error->set_value(ix, 0.0, 0);
           }
        }
      }

      for (size_t j = 0; j < nrefvarsu; j++) {
        size_t i=reflection_ulimit_cv[j];
        for (std::vector<int> ix = target_error->new_index(); target_error->index_ok(ix); target_error->incr(ix) ) {
           cvm::real cv_value=target_error->bin_to_value_scalar(ix[i], i);
           if (cv_value > reflection_ulimit[j]){
             target_error->set_value(ix, 0.0, 0);
           }
        }
      }

      target_error->multiply_constant(intval);

      if (!use_grids) { 
        cvm::error("Error: You must enable grids when using ebMeta with error"
                   "on the target distribution ", INPUT_ERROR);
      }
      bin_volume = 1.0;
      for (size_t i = 0; i < num_variables(); i++) {
          bin_volume*=variables(i)->width;
      }
      target_error->multiply_constant(bin_volume);

      eff_error_points=0;
      for (size_t i = 0; i < npoints; i++) {
         if ( target_error->array_value(i)>0 || target_dist->array_value(i)>0 ) {
           eff_error_points++;
         } 
      }
      if (which_error_point.size()==0) {
        which_error_point.resize(eff_error_points);
      }
      eff_error_points=0;
      for (size_t i = 0; i < npoints; i++) {
         if ( target_error->array_value(i)>0 || target_dist->array_value(i)>0 ) {
           which_error_point[eff_error_points]=i;
           eff_error_points++;
         }         
      }

      gamma_vec.resize(eff_error_points);
      gamma_vec.assign(eff_error_points, 1.0);
      if (target_prob.size()==0) {
        target_prob.resize(eff_error_points);
      }

      if (target_dist_eff.size()==0) {
        target_dist_eff.resize(eff_error_points);
      }

      if (target_error_orig.size()==0) {
        target_error_orig.resize(eff_error_points);
      }

      if (target_prob_orig.size()==0) {
        target_prob_orig.resize(eff_error_points);
      }

      for (size_t j = 0; j < eff_error_points; j++) {
         size_t i=which_error_point[j];
         target_prob[j]=target_dist->array_value(i)*bin_volume;
         target_dist_eff[j]=target_prob[j];
         target_prob_orig[j]=target_prob[j];
         target_error_orig[j]=target_error->array_value(i);  
      }
    }
    // multiply by effective volume = exp(differential entropy)
    cvm::real volume = std::exp(target_dist->entropy());
    ebmeta_tau0 = gauss_factor*volume;
    target_dist->multiply_constant(volume);
    // get and check minimum positive value
    min_pos_val = target_dist->minimum_pos_value();
    if(min_pos_val<=0){
      cvm::error("Error: Target distribution of ebMeta has negative "
                 "or zero minimum positive value!.\n", INPUT_ERROR);
    }
    get_keyval(conf, "ebMetaMaxScaleF", ebmeta_max_scale_f, 1.0/min_pos_val);
    get_keyval(conf, "ebMetaUpdateScaleF", ebmeta_factp, 1.0);
    if (ebmeta_factp/ebmeta_tau0>=1.0) {
      cvm::error("Error: Overall scale factor for target distribution update" 
                 "(ebMetaUpdateScaleF/ebmeta_tau0) is larger than 1.0 please" 
                 "reduce ebMetaUpdateScaleF !.\n", INPUT_ERROR);       
    }
    get_keyval(conf, "ebMetaNormConstToll", ebmeta_nconst_toll, 0.000001);
    get_keyval(conf, "ebMetaNormConstMaxSteps", ebmeta_nconst_maxsteps, 1000);
    get_keyval(conf, "ebMetaUpdateTargets", update_targets, false);
    get_keyval(conf, "ebMetaForgetTargets", forget_targets, false);
    ebmeta_ftarget=1;
    if (forget_targets) {
      update_targets=true;
      ebmeta_ftarget=-1; 
      get_keyval(conf, "ebMetaMaxErrorScale",ebmeta_maxerror_s, 4.0);
    }
    if (update_targets) {
      get_keyval(conf, "ebMetaUpdateErrorScale", update_error_s, 0.1);
      get_keyval(conf, "ebMetaUpdateProbScale", update_prob_s, 0.01);
      get_keyval(conf, "ebMetaMinErrorScale",ebmeta_minerror_s, 0.25);
      ebmeta_minerror=ebmeta_minerror_s*bin_volume/(volume); 
    }    
  }

  return COLVARS_OK;
}

int colvarbias_meta::init_inversion_params(std::string const &conf)
{
  bool use_inversion;
  get_keyval(conf, "useHillsInversion", use_inversion, false);
  ninvvarsl=0;
  ninvvarsu=0;
  if (use_inversion) {
    inversion_type = it_monod;
    std::string inversion_type_str;
    get_keyval(conf, "inversionType", inversion_type_str, to_lower_cppstr(std::string("monoDimensional")));
    inversion_type_str = to_lower_cppstr(inversion_type_str);
    if (inversion_type_str == to_lower_cppstr(std::string("monoDimensional"))) {
      inversion_type = it_monod;
    } else if (inversion_type_str == to_lower_cppstr(std::string("multiDimensional"))) {
      inversion_type = it_multid;
    }

    get_keyval(conf, "inversionMaxHillsWeight", inv_max_ww, 5.0);
    get_keyval(conf, "inversionLowLimitNCVs", ninvvarsl, num_variables());
    get_keyval(conf, "inversionUpLimitNCVs", ninvvarsu, num_variables());
    if (inversion_llimit_cv.size()==0) {
      inversion_llimit_cv.resize(ninvvarsl);
      for (size_t i = 0; i < ninvvarsl; i++) {
         inversion_llimit_cv[i]=i;
      }
    }
    if (inversion_ulimit_cv.size()==0) {
      inversion_ulimit_cv.resize(ninvvarsu);
      for (size_t i = 0; i < ninvvarsu; i++) {
         inversion_ulimit_cv[i]=i;
      }
    }
    if(ninvvarsl>0) {
      if (get_keyval(conf, "inversionLowLimitUseCVs", inversion_llimit_cv, inversion_llimit_cv)) {
        if (inversion_llimit.size()==0) {
          inversion_llimit.resize(ninvvarsl);
        }
      } else {
        cvm::log("Using all variables for lower limits of inversion \n");
      }
      if (get_keyval(conf, "inversionLowLimit", inversion_llimit, inversion_llimit)) {

        if (inversion_intl.size()==0) {
          inversion_intl.resize(ninvvarsl);
        }
        if (!get_keyval(conf, "inversionLowLimitRange", inversion_intl, inversion_intl)) {
          inversion_intl.assign(ninvvarsl,6.0);
        }

        if (inversion_ref_intl.size()==0) {
          inversion_ref_intl.resize(ninvvarsl);
        }
        if (!get_keyval(conf, "inversionRefLLimitRange", inversion_ref_intl, inversion_ref_intl)) {
          inversion_ref_intl.assign(ninvvarsl,1.6);
        }

        for (size_t i = 0; i < ninvvarsl; i++) {
           if (use_grids) {
             size_t ii=inversion_llimit_cv[i];
             cvm:: real sigma=0.5*variables(ii)->width*hill_width;
             cvm:: real bound=variables(ii)->lower_boundary;
             cvm:: real inv_r=inversion_llimit[i]-inversion_intl[i]*sigma;
             if (inv_r < bound) {
               cvm::error("Error: When using grids, lower boundary for CV"+cvm::to_str(ii)+" must be smaller than"+cvm::to_str(inv_r)+".\n", INPUT_ERROR);
             }
           }
           cvm::log("Inversion condition is applied on a lower limit for CV "+cvm::to_str(inversion_llimit_cv[i])+".\n");
           cvm::log("Inversion condition lower limit for this CV is "+cvm::to_str(inversion_llimit[i])+".\n");
           cvm::log("Inversion condition lower range for this CV is "+cvm::to_str(inversion_intl[i])+".\n");
           cvm::log("Inversion condition lower reflection range for this CV is "+cvm::to_str(inversion_ref_intl[i])+".\n");
        }      
      } else {
        cvm::error("Error: Lower limits for inversion not provided.\n", INPUT_ERROR);
        return INPUT_ERROR;
      }
    }
    
    if(ninvvarsu>0) {
      if (get_keyval(conf, "inversionUpLimitUseCVs", inversion_ulimit_cv, inversion_ulimit_cv)) {
        if (inversion_ulimit.size()==0) {
          inversion_ulimit.resize(ninvvarsu);
        }
      } else {
        cvm::log("Using all variables for upper limits of inversion \n");
      }

      if (get_keyval(conf, "inversionUpLimit", inversion_ulimit, inversion_ulimit)) {
        if (inversion_intu.size()==0) {
          inversion_intu.resize(ninvvarsu);
        }
        if (!get_keyval(conf, "inversionUpLimitRange", inversion_intu, inversion_intu)) {
          inversion_intu.assign(ninvvarsu,6.0);
        } 

        if (inversion_ref_intu.size()==0) {
          inversion_ref_intu.resize(ninvvarsu);
        }
        if (!get_keyval(conf, "inversionRefULimitRange", inversion_ref_intu, inversion_ref_intu)) {
          inversion_ref_intu.assign(ninvvarsu,1.6);
        }


        for (size_t i = 0; i < ninvvarsu; i++) {
           if (use_grids) {
             size_t ii=inversion_ulimit_cv[i];
             cvm:: real sigma=0.5*variables(ii)->width*hill_width;
             cvm:: real bound=variables(ii)->upper_boundary;
             cvm:: real inv_r=inversion_ulimit[i]+inversion_intu[i]*sigma;
             if (inv_r > bound) {
               cvm::error("Error: When using grids, upper boundary for CV"+cvm::to_str(ii)+" must be larger than"+cvm::to_str(inv_r)+".\n", INPUT_ERROR);
             }
           }

           cvm::log("Inversion condition is applied on an upper limit for CV "+cvm::to_str(inversion_ulimit_cv[i])+".\n");
           cvm::log("Inversion condition upper limit for this CV is "+cvm::to_str(inversion_ulimit[i])+".\n");
           cvm::log("Inversion condition upper range for this CV is "+cvm::to_str(inversion_intu[i])+".\n");
           cvm::log("Inversion condition upper reflection range for this CV is "+cvm::to_str(inversion_ref_intu[i])+".\n");
        }
      } else {
        cvm::error("Error: Upper limits for inversion not provided.\n", INPUT_ERROR); 
        return INPUT_ERROR;
      }
    }
  }
  // use inversion only for scalar variables

  for (size_t i = 0; i < ninvvarsl; i++) {
     if (inversion_llimit_cv[i]>=num_variables() || inversion_llimit_cv[i]<0) {
       cvm::error("Error: CV number is negative or >= num_variables  \n", INPUT_ERROR);
       return INPUT_ERROR;
     }
     int j=inversion_llimit_cv[i];
     if (variables(j)->value().type()!=colvarvalue::type_scalar) {
       cvm::error("Error: Hills inversion can be used only with scalar variables.\n", INPUT_ERROR);
       return INPUT_ERROR;
     }
  } 

  for (size_t i = 0; i < ninvvarsu; i++) {
     if (inversion_ulimit_cv[i]>=num_variables() || inversion_ulimit_cv[i]<0) {
       cvm::error("Error: CV number is negative or >= num_variables  \n", INPUT_ERROR);
       return INPUT_ERROR;
     }
     int j=inversion_ulimit_cv[i];
     if (variables(j)->value().type()!=colvarvalue::type_scalar) {
       cvm::error("Error: Hills inversion can be used only with scalar variables.\n", INPUT_ERROR);
       return INPUT_ERROR;
     }
  }

  // mono vs multimensional inversion

  switch (inversion_type) {
  case it_monod:
    cvm::log("Using monodimensional inversion \n");
    break;
  case it_multid:
    cvm::log("Using multidimensional inversion \n");
    int sum=1;
    int nstates;
    if (inv_state_ll.size()==0) {
      inv_state_ll.resize(ninvvarsl,std::vector<int>(1));
    }
    inv_state_ll[0][0]=1;
    for (size_t j = 1; j < ninvvarsl; j++) {
      sum*=10;
      nstates=0;
      for (size_t jj = 0; jj < j; jj++) {
            nstates+=inv_state_ll[j].size();
      }
      nstates++;
      inv_state_ll[j].resize(nstates);
      inv_state_ll[j][0]=sum;
      int count=0;
      for (size_t jj = 0; jj < j; jj++) {
         for (size_t ii = 0; ii < inv_state_ll[jj].size(); ii++) {
            count++;
            inv_state_ll[j][count]=inv_state_ll[j][0]+inv_state_ll[jj][ii];
         }
      }
    }

    sum=1;
    if (inv_state_ul.size()==0) {
      inv_state_ul.resize(ninvvarsu,std::vector<int>(1));
    }
    inv_state_ul[0][0]=1;
    for (size_t j = 1; j < ninvvarsu; j++) {
      sum*=10;
      nstates=0;
      for (size_t jj = 0; jj < j; jj++) {
            nstates+=inv_state_ul[j].size();
      }
      nstates++;
      inv_state_ul[j].resize(nstates);
      inv_state_ul[j][0]=sum;
      int count=0;
      for (size_t jj = 0; jj < j; jj++) {
         for (size_t ii = 0; ii < inv_state_ul[jj].size(); ii++) {
            count++;
            inv_state_ul[j][count]=inv_state_ul[j][0]+inv_state_ul[jj][ii];
         }
      }
    }
    break;
  }


  return COLVARS_OK;
}

int colvarbias_meta::init_reflection_params(std::string const &conf)
{
  bool use_reflection;
  nrefvarsl=0;
  nrefvarsu=0;
  get_keyval(conf, "useHillsReflection", use_reflection, false); 
  if (use_reflection) {

    reflection_type = rt_monod;
    std::string reflection_type_str;
    get_keyval(conf, "reflectionType", reflection_type_str, to_lower_cppstr(std::string("monoDimensional")));
    reflection_type_str = to_lower_cppstr(reflection_type_str);
    if (reflection_type_str == to_lower_cppstr(std::string("monoDimensional"))) {
      reflection_type = rt_monod;
    } else if (reflection_type_str == to_lower_cppstr(std::string("multiDimensional"))) {
      reflection_type = rt_multid;
    }

    get_keyval(conf, "reflectionLowLimitNCVs", nrefvarsl, num_variables());
    get_keyval(conf, "reflectionUpLimitNCVs", nrefvarsu, num_variables());
    if (reflection_llimit_cv.size()==0) {
      reflection_llimit_cv.resize(nrefvarsl);
      for (size_t i = 0; i < nrefvarsl; i++) {
         reflection_llimit_cv[i]=i;
      }
    }
    if (reflection_ulimit_cv.size()==0) {
      reflection_ulimit_cv.resize(nrefvarsu);
      for (size_t i = 0; i < nrefvarsu; i++) {
         reflection_ulimit_cv[i]=i;
      }
    }
    if(nrefvarsl>0) {
      if (get_keyval(conf, "reflectionLowLimitUseCVs", reflection_llimit_cv, reflection_llimit_cv)) {
        if (reflection_llimit.size()==0) {
          reflection_llimit.resize(nrefvarsl);
        }
      } else {
        cvm::log("Using all variables for lower limits of reflection \n");
      } 
      if (get_keyval(conf, "reflectionLowLimit", reflection_llimit, reflection_llimit)) {
        if (reflection_intl.size()==0) {
          reflection_intl.resize(nrefvarsl);
        }
        if (!get_keyval(conf, "reflectionLowLimitRange", reflection_intl, reflection_intl)) {
          reflection_intl.assign(nrefvarsl,6.0);
        }
        for (size_t i = 0; i < nrefvarsl; i++) {
           cvm::log("Reflection condition is applied on a lower limit for CV "+cvm::to_str(reflection_llimit_cv[i])+".\n");
           cvm::log("Reflection condition lower limit for this CV is "+cvm::to_str(reflection_llimit[i])+".\n");
           cvm::log("Reflection lower limit range for this CV is "+cvm::to_str(reflection_intl[i])+".\n");
        }
      } else {
        cvm::error("Error: Lower limits for reflection not provided.\n", INPUT_ERROR);
        return INPUT_ERROR;
      }
    } 
    
    if(nrefvarsu>0) {
      if (get_keyval(conf, "reflectionUpLimitUseCVs", reflection_ulimit_cv, reflection_ulimit_cv)) {
        if (reflection_ulimit.size()==0) {
          reflection_ulimit.resize(nrefvarsu);
        }
      } else {
        cvm::log("Using all variables for upper limits of reflection \n");
      }
  
      if (get_keyval(conf, "reflectionUpLimit", reflection_ulimit, reflection_ulimit)) {
        if (reflection_intu.size()==0) {
          reflection_intu.resize(nrefvarsu);
        }
        if (!get_keyval(conf, "reflectionUpLimitRange", reflection_intu, reflection_intu)){
          reflection_intu.assign(nrefvarsu,6.0);   
        }
        for (size_t i = 0; i < nrefvarsu; i++) {
           cvm::log("Reflection condition is applied on an upper limit for CV "+cvm::to_str(reflection_ulimit_cv[i])+".\n");
           cvm::log("Reflection condition upper limit for this CV is "+cvm::to_str(reflection_ulimit[i])+".\n");
           cvm::log("Reflection upper limit range for this CV is "+cvm::to_str(reflection_intu[i])+".\n");
        }
      } else {
        cvm::error("Error: Upper limits for reflection not provided.\n", INPUT_ERROR);
        return INPUT_ERROR; 
      }
    }
  }
  // use reflection only with scalar variables

  for (size_t i = 0; i < nrefvarsl; i++) {
     if (reflection_llimit_cv[i]>=num_variables() || reflection_llimit_cv[i]<0) {
       cvm::error("Error: CV number is negative or >= num_variables  \n", INPUT_ERROR);
       return INPUT_ERROR;
     }
     int j=reflection_llimit_cv[i];
     if (variables(j)->value().type()!=colvarvalue::type_scalar) {
       cvm::error("Error: Hills reflection can be used only with scalar variables.\n", INPUT_ERROR);
       return INPUT_ERROR;
     }
  }

  for (size_t i = 0; i < nrefvarsu; i++) {
     if (reflection_ulimit_cv[i]>=num_variables() || reflection_ulimit_cv[i]<0) {
       cvm::error("Error: CV number is negative or >= num_variables  \n", INPUT_ERROR);
       return INPUT_ERROR;
     }
     int j=reflection_ulimit_cv[i];
     if (variables(j)->value().type()!=colvarvalue::type_scalar) {
       cvm::error("Error: Hills reflection can be used only with scalar variables.\n", INPUT_ERROR);
       return INPUT_ERROR;
     }
  }

  // mono vs multimensional reflection

  switch (reflection_type) {
  case rt_monod:
    cvm::log("Using monodimensional reflection \n");
    break;
  case rt_multid:
    // generate reflection states
    cvm::log("Using multidimensional reflection \n");
    int sum=1;
    int nstates;
    if (ref_state_ll.size()==0) {
      ref_state_ll.resize(nrefvarsl,std::vector<int>(1));
    } 
    ref_state_ll[0][0]=1;
    for (size_t j = 1; j < nrefvarsl; j++) {
      sum*=10;
      nstates=0;
      for (size_t jj = 0; jj < j; jj++) {
            nstates+=ref_state_ll[j].size();
      }
      nstates++;
      ref_state_ll[j].resize(nstates);
      ref_state_ll[j][0]=sum;
      int count=0;
      for (size_t jj = 0; jj < j; jj++) {
         for (size_t ii = 0; ii < ref_state_ll[jj].size(); ii++) {
            count++;
            ref_state_ll[j][count]=ref_state_ll[j][0]+ref_state_ll[jj][ii];
         }
      }
    }

    sum=1;
    if (ref_state_ul.size()==0) {
      ref_state_ul.resize(nrefvarsu,std::vector<int>(1));
    }
    ref_state_ul[0][0]=1;
    for (size_t j = 1; j < nrefvarsu; j++) {
      sum*=10;
      nstates=0;
      for (size_t jj = 0; jj < j; jj++) {
            nstates+=ref_state_ul[j].size();
      }
      nstates++;
      ref_state_ul[j].resize(nstates);
      ref_state_ul[j][0]=sum;
      int count=0;
      for (size_t jj = 0; jj < j; jj++) {
         for (size_t ii = 0; ii < ref_state_ul[jj].size(); ii++) {
            count++;
            ref_state_ul[j][count]=ref_state_ul[j][0]+ref_state_ul[jj][ii];
         }
      }
    }

    break;
  }

  return COLVARS_OK;
}

int colvarbias_meta::init_interval_params(std::string const &conf)
{
  bool use_interval;
  use_interval=false;
  nintvarsl=0;
  nintvarsu=0;
  std::vector<int> interval_llimit_cv;
  std::vector<int> interval_ulimit_cv;
  if (get_keyval(conf, "useHillsInterval", use_interval, use_interval)) {
    if (use_interval) {
      get_keyval(conf, "intervalLowLimitNCVs", nintvarsl, num_variables());
      get_keyval(conf, "intervalUpLimitNCVs", nintvarsu, num_variables());
      interval_llimit_cv.resize(nintvarsl); 
      for (size_t i = 0; i < nintvarsl; i++) {
         interval_llimit_cv[i]=i;
      }
      interval_ulimit_cv.resize(nintvarsu);
      for (size_t i = 0; i < nintvarsu; i++) {
         interval_ulimit_cv[i]=i;
      }
      if(nintvarsl>0) {
        if (get_keyval(conf, "intervalLowLimitUseCVs", interval_llimit_cv, interval_llimit_cv)) {
          if (interval_llimit.size()==0) {
            interval_llimit.resize(nintvarsl);
          }
        } else {
          cvm::log("Using all variables for lower limits of interval \n");
        }
        if (get_keyval(conf, "intervalLowLimit", interval_llimit, interval_llimit)) {
          for (size_t i = 0; i < nintvarsl; i++) {
             cvm::log("Hills forces will be removed beyond a lower limit for CV "+cvm::to_str(interval_llimit_cv[i])+".\n");
             cvm::log("Interval condition lower limit for this CV is "+cvm::to_str(interval_llimit[i])+".\n");
          }     
        } else {
          cvm::error("Error: Lower limits for interval not provided.\n", INPUT_ERROR);
          return INPUT_ERROR;
        }
      }
   
      if(nintvarsu>0) {
        if (get_keyval(conf, "intervalUpLimitUseCVs", interval_ulimit_cv, interval_ulimit_cv)) {
          if (interval_ulimit.size()==0) {
            interval_ulimit.resize(nintvarsu);
          }
        } else {
          cvm::log("Using all variables for upper limits of interval \n");
        }
   
        if (get_keyval(conf, "intervalUpLimit", interval_ulimit, interval_ulimit)) {
          for (size_t i = 0; i < nintvarsu; i++) {
             cvm::log("Hills forces will be removed beyond an upper limit for CV "+cvm::to_str(interval_ulimit_cv[i])+".\n");
             cvm::log("Interval condition upper limit for this CV is "+cvm::to_str(interval_ulimit[i])+".\n");
          }
        } else {
          cvm::error("Error: Upper limits for interval not provided.\n", INPUT_ERROR);
          return INPUT_ERROR;
        }
      }
    }
  } else {
    if (nrefvarsl>0 || nrefvarsu>0) {
      cvm::log("Reflection active: Using by default reflection variables and limits for interval \n");
      nintvarsl=nrefvarsl;
      nintvarsu=nrefvarsu;
      interval_llimit_cv.resize(nintvarsl);
      if (interval_llimit.size()==0) {
        interval_llimit.resize(nintvarsl);
      }
      for (size_t i = 0; i < nintvarsl; i++) {
         interval_llimit_cv[i]=reflection_llimit_cv[i];
         interval_llimit[i]=reflection_llimit[i];
      }
      interval_ulimit_cv.resize(nintvarsu);
      if (interval_ulimit.size()==0) {
        interval_ulimit.resize(nintvarsu);
      }
      for (size_t i = 0; i < nintvarsu; i++) {
         interval_ulimit_cv[i]=reflection_ulimit_cv[i];
         interval_ulimit[i]=reflection_ulimit[i];
      }
    }     
  }

  if (which_int_llimit_cv.size()==0) {
    which_int_llimit_cv.resize(num_variables());
  }
  for (size_t i = 0; i < num_variables(); i++) { 
     which_int_llimit_cv[i]=-1;
  }
  for (size_t i = 0; i < nintvarsl; i++) { 
     int j=interval_llimit_cv[i];
     which_int_llimit_cv[j]=i;
  }

  if (which_int_ulimit_cv.size()==0) { 
    which_int_ulimit_cv.resize(num_variables());
  }
  for (size_t i = 0; i < num_variables(); i++) { 
     which_int_ulimit_cv[i]=-1;
  }
  for (size_t i = 0; i < nintvarsu; i++) { 
     int j=interval_ulimit_cv[i];
     which_int_ulimit_cv[j]=i;
  }

  // use interval only with scalar variables

  for (size_t i = 0; i < nintvarsl; i++) {
     if (interval_llimit_cv[i]>=num_variables() || interval_llimit_cv[i]<0) {
       cvm::error("Error: CV number is negative or >= num_variables  \n", INPUT_ERROR);
       return INPUT_ERROR;
     }
     int j=interval_llimit_cv[i];
     if (variables(j)->value().type()!=colvarvalue::type_scalar) {
       cvm::error("Error: Hills interval can be used only with scalar variables.\n", INPUT_ERROR);
       return INPUT_ERROR;
     }
  }

  for (size_t i = 0; i < nintvarsu; i++) {
     if (interval_ulimit_cv[i]>=num_variables() || interval_ulimit_cv[i]<0) {
       cvm::error("Error: CV number is negative or >= num_variables  \n", INPUT_ERROR);
       return INPUT_ERROR;
     }
     int j=interval_ulimit_cv[i];
     if (variables(j)->value().type()!=colvarvalue::type_scalar) {
       cvm::error("Error: Hills interval can be used only with scalar variables.\n", INPUT_ERROR);
       return INPUT_ERROR;
     }
  }

  return COLVARS_OK;
}

colvarbias_meta::~colvarbias_meta()
{
  if (hills_energy) {
    delete hills_energy;
    hills_energy = NULL;
  }

  if (hills_energy_gradients) {
    delete hills_energy_gradients;
    hills_energy_gradients = NULL;
  }

  if (replica_hills_os) {
    cvm::proxy->close_output_stream(replica_hills_file);
    replica_hills_os = NULL;
  }

  if (hills_traj_os) {
    cvm::proxy->close_output_stream(hills_traj_file_name());
    hills_traj_os = NULL;
  }

  if(target_dist) {
    delete target_dist;
    target_dist = NULL;
  }

  if(target_error) {
    delete target_error;
    target_error = NULL;
  }

  if (ebmeta_out && ebmetaoutfile.is_open())
    ebmetaoutfile.close();

}



// **********************************************************************
// Hill management member functions
// **********************************************************************

std::list<colvarbias_meta::hill>::const_iterator
colvarbias_meta::create_hill(colvarbias_meta::hill const &h)
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
    if (min_dist < (3.0 * std::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(h);
    }
  }

  // output to trajectory (if specified)
  if (hills_traj_os) {
    *hills_traj_os << (hills.back()).output_traj();
    cvm::proxy->flush_output_stream(hills_traj_os);
  }

  has_data = true;
  return hills.end();
}

bool colvarbias_meta::check_reflection_limits(bool &ah)
{
  for (size_t i = 0; i < nrefvarsl; i++) {
     int ii=reflection_llimit_cv[i];
     cvm:: real cv_value=variables(ii)->value();
     if (cv_value<reflection_llimit[i]) {
       ah=false;
     }
  }
  for (size_t i = 0; i < nrefvarsu; i++) {
     int ii=reflection_ulimit_cv[i];
     cvm:: real cv_value=variables(ii)->value();
     if (cv_value>reflection_ulimit[i]) {
       ah=false;
     } 
  }
  return ah;  
}

bool colvarbias_meta::check_inversion_limits(bool &ah)
{
  for (size_t i = 0; i < ninvvarsl; i++) {
     int ii=inversion_llimit_cv[i];
     cvm:: real cv_value=variables(ii)->value();
     if (cv_value<inversion_llimit[i]) {
       ah=false;
     }
  }

  for (size_t i = 0; i < ninvvarsu; i++) {
     int ii=inversion_ulimit_cv[i];
     cvm:: real cv_value=variables(ii)->value();
     if (cv_value>inversion_ulimit[i]) {
       ah=false;
     }
  }
  return ah;
}

int colvarbias_meta::reflect_hill_multid(int const &aa,
                                     cvm::real const &h_scale, 
                                     std::vector<std::vector<int>> const &ref_state,
                                     std::vector<int> const &ref_lim_cv,
                                     std::vector<cvm::real> const &ref_lim,
                                     std::vector<cvm::real> const &ref_int)
{
  size_t i = 0;
  std::vector<colvarvalue> curr_cv_values(num_variables());
  for (i = 0; i < num_variables(); i++) {
    curr_cv_values[i].type(variables(i)->value());
  }
  std::vector<cvm::real> h_w(num_variables());
  for (i = 0; i < num_variables(); i++) {
      curr_cv_values[i] = variables(i)->value();
      h_w[i]=variables(i)->width*hill_width;
  }
  for (size_t j = 0; j < aa; j++) {
     int startsum=1;
     for (size_t i = 0; i < j; i++) {
        startsum*=10;
     }
     for (size_t jj = 0; jj < ref_state[j].size(); jj++) {
        bool hill_add=false;
        int getsum=startsum;
        int countstate=0;
        int check_val=ref_state[j][jj];
        for (size_t i = 0; i <= j; i++) {
           int upordown=std::floor(check_val/getsum);
           int state=aa-1-j+countstate;
           countstate++;
           check_val=check_val-getsum;
           getsum=getsum/10;

           size_t ii=ref_lim_cv[state];
           cvm:: real tmps=0.5*h_w[ii];
           colvarvalue tmp=curr_cv_values[ii]; // store original current cv value
           colvarvalue unitary=curr_cv_values[ii];
           unitary.set_to_one();
           cvm:: real tmpd=ref_lim[state]-cvm::real(curr_cv_values[ii]);
           tmpd=std::sqrt(tmpd*tmpd);
           if (tmpd<ref_int[state]*tmps*upordown ) { // do mirror within selected range in case upordown=1
             hill_add=true;
             curr_cv_values[ii]=2.0*ref_lim[state]*unitary-tmp; // reflected cv value
           }
        }
        if (hill_add) {
          std::string h_replica = "";
          switch (comm) {
       
          case single_replica:
       
            create_hill(hill(cvm::step_absolute(), hill_weight*h_scale, curr_cv_values, h_w, h_replica));
       
            break;
       
          case multiple_replicas:
            h_replica=replica_id;
            create_hill(hill(cvm::step_absolute(), hill_weight*h_scale, curr_cv_values, h_w, replica_id));
            if (replica_hills_os) {
              *replica_hills_os << hills.back();
            } else {
              return cvm::error("Error: in metadynamics bias \""+this->name+"\""+
                                ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                                " while writing hills for the other replicas.\n", FILE_ERROR);
            }
            break;
          }
       
          for (size_t i = 0; i < num_variables(); i++) {
             curr_cv_values[i] = variables(i)->value(); // go back to previous values
          }
        }  
     }
  }
  return COLVARS_OK;
}

int colvarbias_meta::reflect_hill_monod(int const &aa,
                                   cvm::real const &h_scale,
                                   cvm::real const &ref_lim,
                                   cvm::real const &ref_int)

{
  size_t i = 0;
  std::vector<colvarvalue> curr_cv_values(num_variables());
  for (i = 0; i < num_variables(); i++) {
    curr_cv_values[i].type(variables(i)->value());
  }
  std::vector<cvm::real> h_w(num_variables());
  for (i = 0; i < num_variables(); i++) {
      curr_cv_values[i] = variables(i)->value();
      h_w[i]=variables(i)->width*hill_width;
  }
  cvm:: real tmps=0.5*h_w[aa];
  colvarvalue tmp=curr_cv_values[aa]; // store original current cv value
  colvarvalue unitary=curr_cv_values[aa];
  unitary.set_to_one();
  cvm:: real tmpd=ref_lim-cvm::real(curr_cv_values[aa]);
  tmpd=std::sqrt(tmpd*tmpd);
  if (tmpd<ref_int*tmps ) { // do mirror within selected range
    curr_cv_values[aa]=2.0*ref_lim*unitary-tmp; // reflected cv value
    std::string h_replica = "";
    switch (comm) {

    case single_replica:

      create_hill(hill(cvm::step_absolute(), hill_weight*h_scale, curr_cv_values, h_w, h_replica));

      break;

    case multiple_replicas:
      h_replica=replica_id;
      create_hill(hill(cvm::step_absolute(), hill_weight*h_scale, curr_cv_values, h_w, h_replica));
      if (replica_hills_os) {
        *replica_hills_os << hills.back();
      } else {
        return cvm::error("Error: in metadynamics bias \""+this->name+"\""+
                          ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                          " while writing hills for the other replicas.\n", FILE_ERROR);
      }
      break;
    }
    curr_cv_values[aa]=tmp; // go back to previous value
  }
  return COLVARS_OK;
}

int colvarbias_meta::invert_hill_multid(int const &aa,
                                     cvm::real const &h_scale,
                                     std::vector<std::vector<int>> const &inv_state,
                                     std::vector<int> const &inv_lim_cv,
                                     std::vector<cvm::real> const &inv_lim,
                                     std::vector<cvm::real> const &ref_int,
                                     std::vector<cvm::real> const &inv_int)
{
  size_t i = 0;
  std::vector<colvarvalue> curr_cv_values(num_variables());
  std::vector<colvarvalue> ref_cv_values(num_variables());
  std::vector<colvarvalue> curr_bound(num_variables());
  for (i = 0; i < num_variables(); i++) {
    curr_cv_values[i].type(variables(i)->value());
    ref_cv_values[i].type(variables(i)->value());
    curr_bound[i].type(variables(i)->value());
  }
  std::vector<cvm::real> h_w(num_variables());
  for (i = 0; i < num_variables(); i++) {
      curr_cv_values[i] = variables(i)->value();
      ref_cv_values[i] = curr_cv_values[i];
      curr_bound[i] = curr_cv_values[i];
      h_w[i]=variables(i)->width*hill_width;
  }
  cvm:: real smooth_func=0;
  for (size_t j = 0; j < aa; j++) {
     int startsum=1;
     for (size_t i = 0; i < j; i++) {
        startsum*=10;
     }
     for (size_t jj = 0; jj < inv_state[j].size(); jj++) {
        bool hill_add=false;
        bool add_invert=false; 
        cvm:: real hills_ww=hill_weight*h_scale;
        int getsum=startsum;
        int countstate=0;
        int check_val=inv_state[j][jj];
        for (size_t i = 0; i <= j; i++) {
           int upordown=std::floor(check_val/getsum);
           int state=aa-1-j+countstate;
           countstate++;
           check_val=check_val-getsum;
           getsum=getsum/10;
           size_t ii=inv_lim_cv[state];
           cvm:: real tmps=0.5*h_w[ii];
           colvarvalue tmp=ref_cv_values[ii]; // store original current cv value
           colvarvalue unitary=ref_cv_values[ii];
           unitary.set_to_one(); 
           cvm:: real tmpd=inv_lim[state]-cvm::real(ref_cv_values[ii]);
           tmpd=std::sqrt(tmpd*tmpd);
           smooth_func=1/(1+pow(tmpd/(inv_int[state]*ref_int[state]*tmps),10)); 
           if (tmpd<ref_int[state]*tmps*upordown ) { // do mirror within selected range in case upordown=1
             ref_cv_values[ii]=2.0*inv_lim[state]*unitary-tmp; // reflected cv value
             curr_bound[ii]=inv_lim[state]*unitary; //boundary cv value
	   } else if (smooth_func*upordown >= 0.01) {
             add_invert=true;     
             ref_cv_values[ii]=2.0*inv_lim[state]*unitary-tmp; // reflected cv value
             curr_bound[ii]=inv_lim[state]*unitary; //boundary cv value              
           }
        }
        if (add_invert) {
          cvm::real tmpF1 = 0.0; // Bias potential in current position
          if (use_grids) {
            std::vector<int> curr_bin_inv = hills_energy->get_colvars_index(curr_cv_values);
            tmpF1 = hills_energy->value(curr_bin_inv);
          } else {
            calc_hills(new_hills_begin, hills.end(), tmpF1, curr_cv_values);
          }
          cvm::real tmpF2 = 0.0; // Bias potential at the border
          if (use_grids) {
            std::vector<int> curr_bin_inv = hills_energy->get_colvars_index(curr_bound);
            tmpF2 = hills_energy->value(curr_bin_inv);
          } else {
            calc_hills(new_hills_begin, hills.end(), tmpF2, curr_bound);
          }
          cvm::real tmpF3 = 0.0; // Bias potential on symmetric position
          if (use_grids) {
            std::vector<int> curr_bin_inv = hills_energy->get_colvars_index(ref_cv_values);
            tmpF3 = hills_energy->value(curr_bin_inv);
          } else {
            calc_hills(new_hills_begin, hills.end(), tmpF3, ref_cv_values);
          }
          hills_ww = (2.*tmpF2-tmpF1-tmpF3)*smooth_func;
          if (std::sqrt(hills_ww*hills_ww)>inv_max_ww*hill_weight) {
            hills_ww = inv_max_ww*hill_weight*(std::sqrt(hills_ww*hills_ww))/hills_ww;
          }
                    
        }
        if (hill_add) {
          std::string h_replica = "";
          switch (comm) {
   
          case single_replica:
   
            create_hill(hill(cvm::step_absolute(), hills_ww, ref_cv_values, h_w, h_replica));
   
            break;
   
          case multiple_replicas:
            h_replica = replica_id;
            create_hill(hill(cvm::step_absolute(), hills_ww, ref_cv_values, h_w, h_replica));
            if (replica_hills_os) {
              *replica_hills_os << hills.back();
            } else {
              return cvm::error("Error: in metadynamics bias \""+this->name+"\""+
                                ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                                " while writing hills for the other replicas.\n", FILE_ERROR);
            }
            break;
          }
   
          for (size_t i = 0; i < num_variables(); i++) {
             ref_cv_values[i] = variables(i)->value(); // go back to previous values
             curr_bound[i] = ref_cv_values[i];
          }
        }
     }
  }
  return COLVARS_OK;
}


int colvarbias_meta::invert_hill_monod(int const &aa,
                                   cvm::real const &h_scale,
                                   cvm::real const &inv_lim,
                                   cvm::real const &ref_int,
                                   cvm::real const &inv_int)
{
  size_t i = 0;
  std::vector<colvarvalue> curr_cv_values(num_variables());
  for (i = 0; i < num_variables(); i++) {
    curr_cv_values[i].type(variables(i)->value());
  }
  std::vector<cvm::real> h_w(num_variables());
  for (i = 0; i < num_variables(); i++) {
      curr_cv_values[i] = variables(i)->value();
      h_w[i]=variables(i)->width*hill_width;
  }

  cvm:: real tmps=0.5*h_w[aa];
  colvarvalue tmp=curr_cv_values[aa]; // store original current cv value
  colvarvalue unitary=curr_cv_values[aa];
  unitary.set_to_one();
  cvm:: real tmpd=inv_lim-cvm::real(curr_cv_values[aa]);
  tmpd=std::sqrt(tmpd*tmpd);
  cvm:: real smooth_func=1/(1+pow(tmpd/(inv_int*ref_int*tmps),10));
  cvm:: real hills_ww=hill_weight*h_scale;
  bool hill_add=false;
  if ( tmpd < ref_int*tmps ) { // just mirror ..
    hill_add=true;
    curr_cv_values[aa]=2.0*inv_lim*unitary-tmp; // reflected cv value
  } else if (smooth_func >= 0.01) { // do inversion condition     
    hill_add=true;
    cvm::real tmpF1 = 0.0; // Bias potential in current position
    if (use_grids) {
      std::vector<int> curr_bin_inv = hills_energy->get_colvars_index(curr_cv_values);
      tmpF1 = hills_energy->value(curr_bin_inv);
    } else {
      calc_hills(new_hills_begin, hills.end(), tmpF1, curr_cv_values);
    }
    curr_cv_values[aa]=inv_lim*unitary;
    cvm::real tmpF2 = 0.0; // Bias potential at the border
    if (use_grids) {
      std::vector<int> curr_bin_inv = hills_energy->get_colvars_index(curr_cv_values);
      tmpF2 = hills_energy->value(curr_bin_inv);
    } else {
      calc_hills(new_hills_begin, hills.end(), tmpF2, curr_cv_values);
    }
    curr_cv_values[aa]=2.0*inv_lim*unitary-tmp; // reflected cv value
    cvm::real tmpF3 = 0.0; // Bias potential on symmetric position
    if (use_grids) {
      std::vector<int> curr_bin_inv = hills_energy->get_colvars_index(curr_cv_values);
      tmpF3 = hills_energy->value(curr_bin_inv);
    } else {
      calc_hills(new_hills_begin, hills.end(), tmpF3, curr_cv_values);
    }
    hills_ww = (2.*tmpF2-tmpF1-tmpF3)*smooth_func;
    if (std::sqrt(hills_ww*hills_ww)>inv_max_ww*hill_weight) {
      hills_ww = inv_max_ww*hill_weight*(std::sqrt(hills_ww*hills_ww))/hills_ww;
    }
  } 
  if (hill_add) {
    std::string h_replica = ""; 
    switch (comm) {
     
    case single_replica:
   
      create_hill(hill(cvm::step_absolute(), hills_ww, curr_cv_values, h_w, h_replica));
   
      break;
   
    case multiple_replicas:
      h_replica = replica_id;
      create_hill(hill(cvm::step_absolute(), hills_ww, curr_cv_values, h_w, h_replica));
      if (replica_hills_os) {
        *replica_hills_os << hills.back();
      } else {
        return cvm::error("Error: in metadynamics bias \""+this->name+"\""+
                          ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                          " while writing hills for the other replicas.\n", FILE_ERROR);
      }
      break;
    }
    curr_cv_values[aa]=tmp; // go back to previous value
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

  if (hills_traj_os) {
    // output to the trajectory
    *hills_traj_os << "# DELETED this hill: "
                   << (hills.back()).output_traj()
                   << "\n";
    cvm::proxy->flush_output_stream(hills_traj_os);
  }

  return hills.erase(h);
}


int colvarbias_meta::update()
{
  int error_code = COLVARS_OK;

  // update base class
  error_code |= colvarbias::update();

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

  error_code |= calc_energy();
  error_code |= calc_forces();

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
        (3 * (size_t) std::floor(hill_width)) + 1;

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

        if (!variables(i)->hard_lower_boundary)
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

        if (!variables(i)->hard_upper_boundary)
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
  // add a new hill if the required time interval has passed
  if ((cvm::step_absolute() % new_hill_freq) == 0 &&
      is_enabled(f_cvb_history_dependent)) {

    if (cvm::debug()) {
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
               ": adding a new hill at step "+cvm::to_str(cvm::step_absolute())+".\n");
    }

    cvm::real hills_scale=1.0;
    cvm::real volume=0.0; // volume of ebmeta target_dist (in case ebmeta=true)  

    if (ebmeta) {
      cvm:: real scale_fact = target_dist->value(target_dist->get_colvars_index());
      if (scale_fact<=0) {
        scale_fact=min_pos_val;
      }
      if (scale_fact*ebmeta_max_scale_f<1) {
        scale_fact=1.0/ebmeta_max_scale_f;
      } 
      hills_scale *= 1.0/scale_fact;
      if(cvm::step_absolute() <= long(ebmeta_equil_steps)) {
        cvm::real const hills_lambda =
          (cvm::real(long(ebmeta_equil_steps) - cvm::step_absolute())) /
          (cvm::real(ebmeta_equil_steps));
        hills_scale = hills_lambda + (1-hills_lambda)*hills_scale;
      }
    }

    if (ebmetaerror && cvm::step_absolute() > long(ebmeta_equil_steps)) {

       // Update probability distribution according to the error
       // todo this first loop on indices to calculate
       // average bias potential

       cvm::real ave_bias=0.0;
       cvm::real hills_energy_sum_rvalues;
       cvm::real norm_val=0.0;
       
       for (size_t j = 0; j < eff_error_points ; j++) {
         size_t i = which_error_point[j];
         hills_energy_sum_rvalues = 0.0;
         cvm:: real prob_val=target_dist->array_value(i);
         hills_energy_sum_rvalues = hills_energy->array_value(i);
         ave_bias+=hills_energy_sum_rvalues*prob_val;
         norm_val+=prob_val;
       }
       ave_bias=ave_bias/norm_val;

       // Now update target distribution

       for (size_t kk = 0; kk < ebmeta_nconst_maxsteps ; kk++) { // loop to calculate scale factor ebmeta_nconst
          cvm::real sumup=0;
          cvm::real sumdown=0;
          for (size_t j = 0; j < eff_error_points ; j++) {
            size_t i = which_error_point[j];
            cvm::real prob_old=target_dist->array_value(i);
            hills_energy_sum_rvalues = hills_energy->array_value(i);
            cvm::real lambda_e = -(hills_energy_sum_rvalues-ave_bias)/(cvm::kt());
            cvm::real lambda_f = lambda_e-ebmeta_nconst;
            cvm::real exp_prob_value = target_prob[j];
            cvm::real error_val = target_error->array_value(i);
            cvm::real sigma_sqe = error_val*error_val/gamma_vec[j];
            if ( prob_old>0 ) { // update scale factor ebmeta_nconst
              sumup+=sigma_sqe*lambda_e-exp_prob_value+prob_old;
              sumdown+=sigma_sqe;
            }
            cvm::real dev=sigma_sqe*lambda_f;
            if ((std::sqrt(dev*dev))/error_val>ebmeta_fix_bound) {
              dev=error_val*ebmeta_fix_bound*(std::sqrt(dev*dev))/dev;
            }
            cvm::real prob_value = exp_prob_value-dev;
            target_dist->set_array_value(i,prob_value);
          }

          // now we must project the updated probability distribution
          // into the probability simplex (positive and normalized to 1)
          // using algorithm from from Wang, Perpinan 2003
                   
          target_dist->simplexproj(eff_error_points, which_error_point);

          // update scale factor ebmeta_nconst

          cvm::real new_ebmeta_nconst=sumup/sumdown;
          cvm::real ave_nconst=std::sqrt(ebmeta_nconst*ebmeta_nconst);
          cvm::real nconst_toll=(ebmeta_nconst-new_ebmeta_nconst)/ave_nconst;
          nconst_toll=std::sqrt(nconst_toll*nconst_toll);
          if (nconst_toll<=ebmeta_nconst_toll) break;
          ebmeta_nconst=new_ebmeta_nconst;        
       }
       // update gammav
       for (size_t j = 0; j < eff_error_points ; j++) {
          size_t i = which_error_point[j];
          hills_energy_sum_rvalues = hills_energy->array_value(i);
          cvm::real error_val = target_error->array_value(i);
          cvm::real lambda_e = -(hills_energy_sum_rvalues-ave_bias)/(cvm::kt());
          cvm::real lambda_f = lambda_e-ebmeta_nconst;
          cvm::real gammaold = gamma_vec[j];
          cvm::real gammanew = error_val*std::sqrt(lambda_f*lambda_f);
          cvm::real deltagamma=gammanew-gammaold;
          gamma_vec[j] = gammaold+ebmeta_factp*(deltagamma/ebmeta_tau0);
          if(gamma_vec[j]<=0.0) gamma_vec[j]=1.0;
       } 
       // update target dist
       for (size_t j = 0; j < eff_error_points ; j++) {
          size_t i = which_error_point[j];
          cvm:: real probo=target_dist_eff[j];
          cvm:: real probn=target_dist->array_value(i);
          cvm:: real deltaprob=probn-probo;
          target_dist_eff[j]=probo+ebmeta_factp*(deltaprob/ebmeta_tau0);
          target_dist->set_array_value(i,target_dist_eff[j]);
       }
       // update error and central distribution

       if (update_targets) {
         for (size_t j = 0; j < eff_error_points ; j++) {
            size_t i = which_error_point[j];
            cvm::real deltaprob=target_dist_eff[j]-target_prob[j]; 
            deltaprob=std::sqrt(deltaprob*deltaprob);
            cvm::real oerror=target_error->array_value(i);
            cvm::real derror=ebmeta_fix_bound*deltaprob-oerror;
            cvm::real minerror=ebmeta_minerror;
            if(ebmeta_ftarget*minerror>target_error_orig[j]) minerror=target_error_orig[j]; 
            cvm::real nerror=oerror+update_error_s*ebmeta_factp*(derror/ebmeta_tau0);
            if(nerror<minerror) nerror=minerror;
            if(target_error_orig[j]<=0.0) nerror=0.0;
            if(ebmeta_ftarget*nerror>target_error_orig[j]) nerror=target_error_orig[j];
            target_error->set_array_value(i,nerror);
            deltaprob=target_prob_orig[j]-target_dist_eff[j];
            deltaprob=std::sqrt(deltaprob*deltaprob);
            cvm::real deltaprobexp=target_dist_eff[j]-target_prob[j];
            if(ebmeta_ftarget*deltaprob>target_error_orig[j]) deltaprobexp=target_prob_orig[j]-target_prob[j];
            target_prob[j]=target_prob[j]+update_prob_s*ebmeta_factp*(deltaprobexp/ebmeta_tau0);      
         }
       }

       if (forget_targets) {
         cvm::real ebmeta_maxerror=ebmeta_maxerror_s*ebmeta_minerror;
         for (size_t j = 0; j < eff_error_points ; j++) {
            size_t i = which_error_point[j];
            cvm::real nerror=target_error->array_value(i);
            if(nerror>ebmeta_maxerror) target_error->set_array_value(i,ebmeta_maxerror); 
         }    
       }

       // now normalize reference distribution and multiply by effective volume

       cvm::real intval = 1.0/target_dist->integral(eff_error_points, which_error_point);
       target_dist->multiply_constant(eff_error_points, which_error_point, intval);

       volume = std::exp(target_dist->entropy(eff_error_points, which_error_point));
       ebmeta_tau0 = gauss_factor*volume;
       if (ebmeta_factp/ebmeta_tau0>=1) {
         cvm::error("Error: Overall scale factor for target distribution update"
                    "(ebMetaUpdateScaleF/ebmeta_tau0) is larger than 1.0 please"
                    "reduce ebMetaUpdateScaleF !.\n", INPUT_ERROR);
       }
       target_dist->multiply_constant(eff_error_points, which_error_point, volume);
       min_pos_val = target_dist->minimum_pos_value(eff_error_points, which_error_point);
       if (update_targets) ebmeta_minerror=ebmeta_minerror_s*bin_volume/volume;


       // DEBUG print updated distribution and average deviation from experimental value  

       cvm:: real averdev=0.0;
       norm_val=0.0;
       for (size_t j = 0; j < eff_error_points ; j++) {
          size_t i = which_error_point[j];
          cvm::real prob_val=target_dist_eff[j];
          cvm::real error_val = target_error->array_value(i);
          cvm:: real dev=(prob_val-target_prob[j])/error_val;
          dev=std::sqrt(dev*dev);
          averdev+=dev*prob_val;
          norm_val+=prob_val;
       }
       averdev=averdev/norm_val;
       if (ebmeta_out) { 
         if (ebmetaoutfile.is_open()){
           for (size_t j = 0; j < eff_error_points ; j++) {
             size_t i = which_error_point[j];  
             ebmetaoutfile << "EBMetaDdist: ";
             ebmetaoutfile << " ";
             ebmetaoutfile << cvm::step_absolute();
             ebmetaoutfile << " ";
             ebmetaoutfile << target_dist_eff[j];       
             ebmetaoutfile << " ";
             ebmetaoutfile << target_error->array_value(i);
             ebmetaoutfile << " ";
             ebmetaoutfile << target_prob[j];
             ebmetaoutfile << " ";
             ebmetaoutfile << gamma_vec[j];
             ebmetaoutfile << " " << "\n";        
        
             //printf("EBMetaD_error: %i %f %f %f %f %f %f \n",cvm::step_absolute(),target_dist->array_value(i),averdev,gamma_vec[j],target_error->array_value(i),target_prob[j],ebmeta_nconst);    

           }
           ebmetaoutfile << "EBMetaDdata: ";
           ebmetaoutfile << " ";
           ebmetaoutfile << cvm::step_absolute();
           ebmetaoutfile << " ";         
           ebmetaoutfile << averdev;
           ebmetaoutfile << " ";
           ebmetaoutfile << ebmeta_nconst;
           ebmetaoutfile << " ";
           ebmetaoutfile << min_pos_val;
           ebmetaoutfile << " ";
           ebmetaoutfile << intval;
           ebmetaoutfile << " ";
           ebmetaoutfile << ebmeta_minerror/bin_volume;
           ebmetaoutfile << " ";
           ebmetaoutfile << volume;
           ebmetaoutfile << " "<< "\n" ;
           ebmetaoutfile.flush();
         }
       }

       //printf("MinPosVal: %i %f %f %f %f \n",cvm::step_absolute(),min_pos_val,intval,volume,gauss_factor);

    }

    if (well_tempered) {
      cvm::real hills_energy_sum_here = 0.0;
      if (use_grids) {
        std::vector<int> curr_bin = hills_energy->get_colvars_index();
        hills_energy_sum_here = hills_energy->value(curr_bin);
      } else {
        calc_hills(new_hills_begin, hills.end(), hills_energy_sum_here);
      }
      hills_scale *= std::exp(-1.0*hills_energy_sum_here/(bias_temperature*cvm::boltzmann()));
    }

    if (scale_kernel) {
      if (ebmeta && default_kernel_ebmeta) { 
        kernel_coupling_time=cvm::dt() * new_hill_freq*ebmeta_tau0;    
      }
      switch (kernel_type) {
      case kt_inv_sqrt_time:
        hills_scale*= 1.0/std::sqrt(1.0 + (cvm::dt()*cvm::step_absolute()/kernel_coupling_time));
        break;
      case kt_none:
      case kt_ntot:
        break;
      }
    }

    // Whether add a hill
    bool add_hill=true;

    // Do not add hills beyond inversion or reflection borders
    // as just reflected or inverted hills must be present
    // beyond those boundaries

    // Check reflection borders: if beyond borders do not add hill

    add_hill=check_reflection_limits(add_hill);

    // Check inversion borders: if beyond borders do not add hill
    
    add_hill=check_inversion_limits(add_hill);

    if (add_hill) {
      switch (comm) {
   
      case single_replica:
   
        create_hill(hill(hill_weight*hills_scale, colvars, hill_width));
   
        break;
   
      case multiple_replicas:
        create_hill(hill(hill_weight*hills_scale, colvars, hill_width, replica_id));
        if (replica_hills_os) {
          *replica_hills_os << hills.back();
        } else {
          return cvm::error("Error: in metadynamics bias \""+this->name+"\""+
                            ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                            " while writing hills for the other replicas.\n", FILE_ERROR);
        }
        break;
      }
   
      // add reflected and/or inverted hills if required

      switch (reflection_type) {
      case rt_monod: 
        for (size_t i = 0; i < nrefvarsl; i++) {
           size_t ii=reflection_llimit_cv[i];
           reflect_hill_monod(ii, hills_scale, reflection_llimit[i], reflection_intl[i]);
        }     
     
        for (size_t i = 0; i < nrefvarsu; i++) {
           size_t ii=reflection_ulimit_cv[i];
           reflect_hill_monod(ii, hills_scale, reflection_ulimit[i], reflection_intu[i]);
        }
        break;
      case rt_multid:
        reflect_hill_multid(nrefvarsl, hills_scale, ref_state_ll, reflection_llimit_cv, reflection_llimit, reflection_intl);
        reflect_hill_multid(nrefvarsu, hills_scale, ref_state_ul, reflection_ulimit_cv ,reflection_ulimit, reflection_intu);
        break;
      }  
      
      switch (inversion_type) {
      case it_monod: 
        for (size_t i = 0; i < ninvvarsl; i++) {
           size_t ii=inversion_llimit_cv[i];
           invert_hill_monod(ii, hills_scale, inversion_llimit[i], inversion_ref_intl[i], inversion_intl[i]);  
        }
     
        for (size_t i = 0; i < ninvvarsu; i++) {
           size_t ii=inversion_ulimit_cv[i];
           invert_hill_monod(ii, hills_scale, inversion_ulimit[i], inversion_ref_intu[i], inversion_intu[i]);   
        }
        break;
      case it_multid:
        invert_hill_multid(ninvvarsl, hills_scale, inv_state_ll, inversion_llimit_cv, inversion_llimit, inversion_ref_intl, inversion_intl);
        invert_hill_multid(ninvvarsu, hills_scale, inv_state_ul, inversion_ulimit_cv, inversion_ulimit, inversion_ref_intu, inversion_intu);
        break;
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
                  hills_energy,    hills_energy_gradients);
    new_hills_begin = hills.end();

    // TODO: we may want to condense all into one replicas array,
    // including "this" as the first element
    if (comm == multiple_replicas) {
      for (size_t ir = 0; ir < replicas.size(); ir++) {
        replicas[ir]->project_hills(replicas[ir]->new_hills_begin,
                                    replicas[ir]->hills.end(),
                                    replicas[ir]->hills_energy,
                                    replicas[ir]->hills_energy_gradients);
        replicas[ir]->new_hills_begin = replicas[ir]->hills.end();
      }
    }
  }

  return COLVARS_OK;
}


int colvarbias_meta::calc_energy(std::vector<colvarvalue> const &values)
{
  size_t ir = 0;

  for (ir = 0; ir < replicas.size(); ir++) {
    replicas[ir]->bias_energy = 0.0;
  }

  std::vector<int> const curr_bin = values.size() ?
    hills_energy->get_colvars_index(values) :
    hills_energy->get_colvars_index();

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
    for (ir = 0; ir < replicas.size(); ir++) {
      calc_hills(replicas[ir]->hills_off_grid.begin(),
                 replicas[ir]->hills_off_grid.end(),
                 bias_energy,
                 values);
    }
  }

  // now include the hills that have not been binned yet (starting
  // from new_hills_begin)

  for (ir = 0; ir < replicas.size(); ir++) {
    calc_hills(replicas[ir]->new_hills_begin,
               replicas[ir]->hills.end(),
               bias_energy);
    if (cvm::debug()) {
      cvm::log("Hills energy = "+cvm::to_str(bias_energy)+".\n");
    }
  }

  return COLVARS_OK;
}


int colvarbias_meta::calc_forces(std::vector<colvarvalue> const &values)
{
  size_t ir = 0, ic = 0;
  for (ir = 0; ir < replicas.size(); ir++) {
    for (ic = 0; ic < num_variables(); ic++) {
      replicas[ir]->colvar_forces[ic].reset();
    }
  }

  std::vector<int> const curr_bin = values.size() ?
    hills_energy->get_colvars_index(values) :
    hills_energy->get_colvars_index();

  if (hills_energy->index_ok(curr_bin)) {
    for (ir = 0; ir < replicas.size(); ir++) {
      cvm::real const *f = &(replicas[ir]->hills_energy_gradients->value(curr_bin));
      for (ic = 0; ic < num_variables(); ic++) {
        // the gradients are stored, not the forces
        colvar_forces[ic].real_value += -1.0 * f[ic];
      }
    }
  } else {
    // off the grid: compute analytically only the hills at the grid's edges
    for (ir = 0; ir < replicas.size(); ir++) {
      for (ic = 0; ic < num_variables(); ic++) {
        calc_hills_force(ic,
                         replicas[ir]->hills_off_grid.begin(),
                         replicas[ir]->hills_off_grid.end(),
                         colvar_forces,
                         values);
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
      calc_hills_force(ic,
                       replicas[ir]->new_hills_begin,
                       replicas[ir]->hills.end(),
                       colvar_forces,
                       values);
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
                                 std::vector<colvarvalue> const &colvar_values)
{
  size_t i = 0;
  std::vector<colvarvalue> curr_values(num_variables());
  for (i = 0; i < num_variables(); i++) {
    curr_values[i].type(variables(i)->value());
  }

  if (colvar_values.size()) {
    for (i = 0; i < num_variables(); i++) {
      curr_values[i] = colvar_values[i];
    }
  } else {
    for (i = 0; i < num_variables(); i++) {
      curr_values[i] = variables(i)->value();
    }
  }

  for (hill_iter h = h_first; h != h_last; h++) {

    // compute the gaussian exponent
    cvm::real cv_sqdev = 0.0;
    for (i = 0; i < num_variables(); i++) {
      colvarvalue const &x  = curr_values[i];
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      cv_sqdev += (variables(i)->dist2(x, center)) / (half_width*half_width);
    }

    // compute the gaussian
    if (cv_sqdev > 23.0) {
      // set it to zero if the exponent is more negative than log(1.0E-05)
      h->value(0.0);
    } else {
      h->value(std::exp(-0.5*cv_sqdev));
    }
    energy += h->energy();
  }
}


void colvarbias_meta::calc_hills_force(size_t const &i,
                                       colvarbias_meta::hill_iter      h_first,
                                       colvarbias_meta::hill_iter      h_last,
                                       std::vector<colvarvalue>       &forces,
                                       std::vector<colvarvalue> const &values)
{
  // Retrieve the value of the colvar
  colvarvalue const x(values.size() ? values[i] : variables(i)->value());

  // do the type check only once (all colvarvalues in the hills series
  // were already saved with their types matching those in the
  // colvars)

  hill_iter h;
  switch (variables(i)->value().type()) {

  case colvarvalue::type_scalar:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      // if outside interval boundaries do not add force
      bool add_force=true;
      int ii=which_int_llimit_cv[i];
      if (ii>-1 && x[i]<interval_llimit[ii] ) {
        add_force=false;
      }  
      ii=which_int_ulimit_cv[i];
      if (ii>-1 && x[i]>interval_ulimit[ii] ) {
        add_force=false;
      }
      if (add_force) {
        forces[i].real_value +=
          ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
            (variables(i)->dist2_lgrad(x, center)).real_value );
      }
    }
    break;

  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].rvector_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (variables(i)->dist2_lgrad(x, center)).rvector_value );
    }
    break;

  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].quaternion_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (variables(i)->dist2_lgrad(x, center)).quaternion_value );
    }
    break;

  case colvarvalue::type_vector:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].vector1d_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
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
                                    bool print_progress)
{
  if (cvm::debug())
    cvm::log("Metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
             ": projecting hills.\n");

  // TODO: improve it by looping over a small subgrid instead of the whole grid

  std::vector<colvarvalue> colvar_values(num_variables());
  std::vector<cvm::real> colvar_forces_scalar(num_variables());

  std::vector<int> he_ix = he->new_index();
  std::vector<int> hg_ix = (hg != NULL) ? hg->new_index() : std::vector<int> (0);
  cvm::real hills_energy_here = 0.0;
  std::vector<colvarvalue> hills_forces_here(num_variables(), 0.0);

  size_t count = 0;
  size_t const print_frequency = ((hills.size() >= 1000000) ? 1 : (1000000/(hills.size()+1)));

  if (hg != NULL) {

    // loop over the points of the grid
    for ( ;
          (he->index_ok(he_ix)) && (hg->index_ok(hg_ix));
          count++) {
      size_t i;
      for (i = 0; i < num_variables(); i++) {
        colvar_values[i] = hills_energy->bin_to_value_scalar(he_ix[i], i);
      }

      // loop over the hills and increment the energy grid locally
      hills_energy_here = 0.0;
      calc_hills(h_first, h_last, hills_energy_here, colvar_values);
      he->acc_value(he_ix, hills_energy_here);

      for (i = 0; i < num_variables(); i++) {
        hills_forces_here[i].reset();
        calc_hills_force(i, h_first, h_last, hills_forces_here, colvar_values);
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

    // simpler version, with just the energy

    for ( ; (he->index_ok(he_ix)); ) {

      for (size_t i = 0; i < num_variables(); i++) {
        colvar_values[i] = hills_energy->bin_to_value_scalar(he_ix[i], i);
      }

      hills_energy_here = 0.0;
      calc_hills(h_first, h_last, hills_energy_here, colvar_values);
      he->acc_value(he_ix, hills_energy_here);

      he->incr(he_ix);

      count++;
      if ((count % print_frequency) == 0) {
        if (print_progress) {
          cvm::real const progress = cvm::real(count) / cvm::real(he->number_of_points());
          std::ostringstream os;
          os.setf(std::ios::fixed, std::ios::floatfield);
          os << std::setw(6) << std::setprecision(2)
             << 100.0 * progress
             << "% done.";
          cvm::log(os.str());
        }
      }
    }
  }

  if (print_progress) {
    cvm::log("100.00% done.");
  }

  if (! keep_hills) {
    hills.erase(hills.begin(), hills.end());
  }
}


void colvarbias_meta::recount_hills_off_grid(colvarbias_meta::hill_iter  h_first,
                                             colvarbias_meta::hill_iter  h_last,
                                             colvar_grid_scalar         *he)
{
  hills_off_grid.clear();

  for (hill_iter h = h_first; h != h_last; h++) {
    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries(h->centers, true);
    if (min_dist < (3.0 * std::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(*h);
    }
  }
}



// **********************************************************************
// multiple replicas functions
// **********************************************************************


int colvarbias_meta::replica_share()
{
  // sync with the other replicas (if needed)
  if (comm == multiple_replicas) {
    // reread the replicas registry
    update_replicas_registry();
    // empty the output buffer
    if (replica_hills_os) {
      cvm::proxy->flush_output_stream(replica_hills_os);
    }
    read_replica_files();
  }
  return COLVARS_OK;
}


void colvarbias_meta::update_replicas_registry()
{
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
      cvm::error("Error: failed to open file \""+replicas_registry_file+
                 "\" for reading.\n", FILE_ERROR);
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
        if (cvm::debug())
          cvm::log("Metadynamics bias \""+this->name+"\""+
                   ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                   ": skipping this replica's own record: \""+
                   new_replica+"\", \""+new_replica_file+"\"\n");
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
      }
    }
  } else {
    cvm::fatal_error("Error: cannot read the replicas registry file \""+
                     replicas_registry+"\".\n");
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
      (replicas[ir])->update_status = 0;
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
}


void colvarbias_meta::read_replica_files()
{
  // Note: we start from the 2nd replica.
  for (size_t ir = 1; ir < replicas.size(); ir++) {

    if (! (replicas[ir])->replica_state_file_in_sync) {
      // if a new state file is being read, the hills file is also new
      (replicas[ir])->replica_hills_file_pos = 0;
    }

    // (re)read the state file if necessary
    if ( (! (replicas[ir])->has_data) ||
         (! (replicas[ir])->replica_state_file_in_sync) ) {

      cvm::log("Metadynamics bias \""+this->name+"\""+
               ": reading the state of replica \""+
               (replicas[ir])->replica_id+"\" from file \""+
               (replicas[ir])->replica_state_file+"\".\n");

      std::ifstream is((replicas[ir])->replica_state_file.c_str());
      if (! (replicas[ir])->read_state(is)) {
        cvm::log("Reading from file \""+(replicas[ir])->replica_state_file+
                 "\" failed or incomplete: will try again in "+
                 cvm::to_str(replica_update_freq)+" steps.\n");
      } else {
        // state file has been read successfully
        (replicas[ir])->replica_state_file_in_sync = true;
        (replicas[ir])->update_status = 0;
      }
      is.close();
    }

    // now read the hills added after writing the state file
    if ((replicas[ir])->replica_hills_file.size()) {

      if (cvm::debug())
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": checking for new hills from replica \""+
                 (replicas[ir])->replica_id+"\" in the file \""+
                 (replicas[ir])->replica_hills_file+"\".\n");

      // read hills from the other replicas' files; for each file, resume
      // the position recorded previously

      std::ifstream is((replicas[ir])->replica_hills_file.c_str());
      if (is.is_open()) {

        // try to resume the previous position
        is.seekg((replicas[ir])->replica_hills_file_pos, std::ios::beg);
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
            //           if (cvm::debug())
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ": received a hill from replica \""+
                     (replicas[ir])->replica_id+
                     "\" at step "+
                     cvm::to_str(((replicas[ir])->hills.back()).it)+".\n");
          }
          is.clear();
          // store the position for the next read
          (replicas[ir])->replica_hills_file_pos = is.tellg();
          if (cvm::debug())
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ": stopped reading file \""+(replicas[ir])->replica_hills_file+
                     "\" at position "+
                     cvm::to_str((replicas[ir])->replica_hills_file_pos)+".\n");

          // test whether this is the end of the file
          is.seekg(0, std::ios::end);
          if (is.tellg() > (replicas[ir])->replica_hills_file_pos+1) {
            (replicas[ir])->update_status++;
          } else {
            (replicas[ir])->update_status = 0;
          }
        }

      } else {
        cvm::log("Failed to read the file \""+(replicas[ir])->replica_hills_file+
                 "\": will try again in "+
                 cvm::to_str(replica_update_freq)+" steps.\n");
        (replicas[ir])->update_status++;
        // cvm::fatal_error ("Error: cannot read from file \""+
        //                   (replicas[ir])->replica_hills_file+"\".\n");
      }
      is.close();
    }

    size_t const n_flush = (replica_update_freq/new_hill_freq + 1);
    if ((replicas[ir])->update_status > 3*n_flush) {
      // TODO: suspend the calculation?
      cvm::log("WARNING: in metadynamics bias \""+this->name+"\""+
               " failed to read completely the output of replica \""+
               (replicas[ir])->replica_id+
               "\" after more than "+
               cvm::to_str((replicas[ir])->update_status * replica_update_freq)+
               " steps.  Ensure that it is still running.\n");
    }
  }
}


int colvarbias_meta::set_state_params(std::string const &state_conf)
{
  std::string new_replica = "";
  if (colvarparse::get_keyval(state_conf, "replicaID", new_replica,
                              std::string(""), colvarparse::parse_silent) &&
      (new_replica != this->replica_id)) {
    cvm::error("Error: in the state file, the "
               "\"metadynamics\" block has a different replicaID: different system?\n",
               INPUT_ERROR);
    return INPUT_ERROR;
  }

  return COLVARS_OK;
}


std::istream & colvarbias_meta::read_state_data(std::istream& is)
{
  bool grids_from_restart_file = use_grids;

  if (use_grids) {

    if (expand_grids) {
      // the boundaries of the colvars may have been changed; TODO:
      // this reallocation is only for backward-compatibility, and may
      // be deleted when grid_parameters (i.e. colvargrid's own
      // internal reallocation) has kicked in
      delete hills_energy;
      delete hills_energy_gradients;
      hills_energy = new colvar_grid_scalar(colvars);
      hills_energy_gradients = new colvar_grid_gradient(colvars);
    }

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

    size_t const hills_energy_pos = is.tellg();
    std::string key;
    if (!(is >> key)) {
      if (hills_energy_backup != NULL) {
        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy           = hills_energy_backup;
        hills_energy_gradients = hills_energy_gradients_backup;
      }
      is.clear();
      is.seekg(hills_energy_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    } else if (!(key == std::string("hills_energy")) ||
               !(hills_energy->read_restart(is))) {
      is.clear();
      is.seekg(hills_energy_pos, std::ios::beg);
      grids_from_restart_file = false;
      if (!rebin_grids) {
        if (hills_energy_backup == NULL)
          cvm::fatal_error("Error: couldn't read the free energy grid for metadynamics bias \""+
                           this->name+"\""+
                           ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                           "; if useGrids was off when the state file was written, "
                           "enable rebinGrids now to regenerate the grids.\n");
        else {
          if (comm == single_replica)
            cvm::log("Error: couldn't read the free energy grid for metadynamics bias \""+
                     this->name+"\".\n");
          delete hills_energy;
          delete hills_energy_gradients;
          hills_energy           = hills_energy_backup;
          hills_energy_gradients = hills_energy_gradients_backup;
          is.setstate(std::ios::failbit);
          return is;
        }
      }
    }

    size_t const hills_energy_gradients_pos = is.tellg();
    if (!(is >> key)) {
      if (hills_energy_backup != NULL)  {
        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy           = hills_energy_backup;
        hills_energy_gradients = hills_energy_gradients_backup;
      }
      is.clear();
      is.seekg(hills_energy_gradients_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    } else if (!(key == std::string("hills_energy_gradients")) ||
               !(hills_energy_gradients->read_restart(is))) {
      is.clear();
      is.seekg(hills_energy_gradients_pos, std::ios::beg);
      grids_from_restart_file = false;
      if (!rebin_grids) {
        if (hills_energy_backup == NULL)
          cvm::fatal_error("Error: couldn't read the free energy gradients grid for metadynamics bias \""+
                           this->name+"\""+
                           ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                           "; if useGrids was off when the state file was written, "
                           "enable rebinGrids now to regenerate the grids.\n");
        else {
          if (comm == single_replica)
            cvm::log("Error: couldn't read the free energy gradients grid for metadynamics bias \""+
                     this->name+"\".\n");
          delete hills_energy;
          delete hills_energy_gradients;
          hills_energy           = hills_energy_backup;
          hills_energy_gradients = hills_energy_gradients_backup;
          is.setstate(std::ios::failbit);
          return is;
        }
      }
    }

    if (cvm::debug())
      cvm::log("Successfully read new grids for bias \""+
               this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+"\n");

    if (hills_energy_backup != NULL) {
      // now that we have successfully updated the grids, delete the
      // backup copies
      if (cvm::debug())
        cvm::log("Deallocating the older grids.\n");

      delete hills_energy_backup;
      delete hills_energy_gradients_backup;
    }
  }

  bool const existing_hills = !hills.empty();
  size_t const old_hills_size = hills.size();
  hill_iter old_hills_end = hills.end();
  hill_iter old_hills_off_grid_end = hills_off_grid.end();

  // read the hills explicitly written (if there are any)
  while (read_hill(is)) {
    if (cvm::debug())
      cvm::log("Read a previously saved hill under the "
               "metadynamics bias \""+
               this->name+"\", created at step "+
               cvm::to_str((hills.back()).it)+".\n");
  }
  is.clear();
  new_hills_begin = hills.end();
  if (grids_from_restart_file) {
    if (hills.size() > old_hills_size)
      cvm::log("Read "+cvm::to_str(hills.size())+
               " hills in addition to the grids.\n");
  } else {
    if (!hills.empty())
      cvm::log("Read "+cvm::to_str(hills.size())+" hills.\n");
  }

  if (rebin_grids) {

    // allocate new grids (based on the new boundaries and widths just
    // read from the configuration file), and project onto them the
    // grids just read from the restart file

    colvar_grid_scalar   *new_hills_energy =
      new colvar_grid_scalar(colvars);
    colvar_grid_gradient *new_hills_energy_gradients =
      new colvar_grid_gradient(colvars);

    if (!grids_from_restart_file || (keep_hills && !hills.empty())) {
      // if there are hills, recompute the new grids from them
      cvm::log("Rebinning the energy and forces grids from "+
               cvm::to_str(hills.size())+" hills (this may take a while)...\n");
      project_hills(hills.begin(), hills.end(),
                    new_hills_energy, new_hills_energy_gradients, true);
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

  if (use_grids) {
    if (!hills_off_grid.empty()) {
      cvm::log(cvm::to_str(hills_off_grid.size())+" hills are near the "
               "grid boundaries: they will be computed analytically "
               "and saved to the state files.\n");
    }
  }

  if (ebmeta && ebmetaerror) {
    read_ebmetaerror(is);
  }

  if (cvm::debug())
    cvm::log("colvarbias_meta::read_restart() done\n");

  if (existing_hills) {
    hills.erase(hills.begin(), old_hills_end);
    hills_off_grid.erase(hills_off_grid.begin(), old_hills_off_grid_end);
  }

  has_data = true;

  if (comm != single_replica) {
    read_replica_files();
  }

  return is;
}

std::istream & colvarbias_meta::read_hill(std::istream &is)
{
  if (!is) return is; // do nothing if failbit is set

  size_t const start_pos = is.tellg();

  std::string data;
  if ( !(is >> read_block("hill", data)) ) {
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  size_t h_it;
  get_keyval(data, "step", h_it, 0, parse_silent);
  if (h_it <= state_file_step) {
    if (cvm::debug())
      cvm::log("Skipping a hill older than the state file for metadynamics bias \""+
               this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+"\n");
    return is;
  }

  cvm::real h_weight;
  get_keyval(data, "weight", h_weight, hill_weight, parse_silent);

  std::vector<colvarvalue> h_centers(num_variables());
  for (size_t i = 0; i < num_variables(); i++) {
    h_centers[i].type(variables(i)->value());
  }
  {
    // it is safer to read colvarvalue objects one at a time;
    // TODO: change this it later
    std::string centers_input;
    key_lookup(data, "centers", &centers_input);
    std::istringstream centers_is(centers_input);
    for (size_t i = 0; i < num_variables(); i++) {
      centers_is >> h_centers[i];
    }
  }

  std::vector<cvm::real> h_widths(num_variables());
  get_keyval(data, "widths", h_widths,
             std::vector<cvm::real>(num_variables(), (std::sqrt(2.0 * PI) / 2.0)),
             parse_silent);

  std::string h_replica = "";
  if (comm != single_replica) {
    get_keyval(data, "replicaID", h_replica, replica_id, parse_silent);
    if (h_replica != replica_id)
      cvm::fatal_error("Error: trying to read a hill created by replica \""+h_replica+
                       "\" for replica \""+replica_id+
                       "\"; did you swap output files?\n");
  }

  hill_iter const hills_end = hills.end();
  hills.push_back(hill(h_it, h_weight, h_centers, h_widths, h_replica));
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
    if (min_dist < (3.0 * std::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(hills.back());
    }
  }

  has_data = true;
  return is;
}

std::istream & colvarbias_meta::read_ebmetaerror(std::istream &is)
{
  if (!is) return is; // do nothing if failbit is set

  size_t const start_pos = is.tellg();

  std::string data;
  if ( !(is >> read_block("ebmeta", data)) ) {
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  if (!get_keyval(data, "ebmetaNconst", ebmeta_nconst))
    cvm::error("Error: ebmetaNconst missing from the restart.\n");
  if (!get_keyval(data, "targetDist", target_dist_eff))
    cvm::error("Error: targetDist missing from the restart.\n");
  for (size_t j = 0; j < eff_error_points; j++) {
     size_t i = which_error_point[j];
     target_dist->set_array_value(i,target_dist_eff[j]);
  }
  cvm::real intval = 1.0/target_dist->integral(eff_error_points, which_error_point);
  target_dist->multiply_constant(eff_error_points, which_error_point, intval);

  cvm::real volume = std::exp(target_dist->entropy(eff_error_points, which_error_point));
  ebmeta_tau0 = gauss_factor*volume;
  if (ebmeta_factp/ebmeta_tau0>=1) {
    cvm::error("Error: Overall scale factor for target distribution update"
               "(ebMetaUpdateScaleF/ebmeta_tau0) is larger than 1.0 please"
               "reduce ebMetaUpdateScaleF !.\n", INPUT_ERROR);
  }
  target_dist->multiply_constant(eff_error_points, which_error_point, volume);
  min_pos_val = target_dist->minimum_pos_value(eff_error_points, which_error_point);
  if (!get_keyval(data, "gammaVec", gamma_vec))
    cvm::error("Error: gammaVec missing from the restart.\n");
  if (update_targets) {
    std::vector<cvm::real> target_error_tmp(eff_error_points);
    if (!get_keyval(data, "targetError", target_error_tmp))
      cvm::error("Error: targetError missing from the restart.\n");
      for (size_t j = 0; j < eff_error_points; j++) {
         size_t i = which_error_point[j];
         target_error->set_array_value(i,target_error_tmp[j]);
      }
    if (!get_keyval(data, "targetProb", target_prob))
      cvm::error("Error: targetProb missing from the restart.\n");
    ebmeta_minerror=ebmeta_minerror_s*bin_volume/volume;
  }
  return is;
}

int colvarbias_meta::setup_output()
{
  output_prefix = cvm::output_prefix();
  if (cvm::num_biases_feature(colvardeps::f_cvb_calc_pmf) > 1) {
    // if this is not the only free energy integrator, append
    // this bias's name, to distinguish it from the output of the other
    // biases producing a .pmf file
    output_prefix += ("."+this->name);
  }

  if (comm == multiple_replicas) {

    // TODO: one may want to specify the path manually for intricated filesystems?
    char *pwd = new char[3001];
    if (GETCWD(pwd, 3000) == NULL)
      cvm::fatal_error("Error: cannot get the path of the current working directory.\n");
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
    if (!replica_hills_os) {
      cvm::proxy->backup_file(replica_hills_file);
      replica_hills_os = cvm::proxy->output_stream(replica_hills_file);
      if (!replica_hills_os) return cvm::get_error();
      replica_hills_os->setf(std::ios::scientific, std::ios::floatfield);
    }

    // write the state file (so that there is always one available)
    write_replica_state_file();
    // schedule to read the state files of the other replicas
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      (replicas[ir])->replica_state_file_in_sync = false;
    }

    // if we're running without grids, use a growing list of "hills" files
    // otherwise, just one state file and one "hills" file as buffer
    std::ostream *list_os =
      cvm::proxy->output_stream(replica_list_file,
                                (use_grids ? std::ios_base::trunc :
                                 std::ios_base::app));
    if (!list_os) {
      return cvm::get_error();
    }
    *list_os << "stateFile " << replica_state_file << "\n";
    *list_os << "hillsFile " << replica_hills_file << "\n";
    cvm::proxy->close_output_stream(replica_list_file);

    // finally, add a new record for this replica to the registry
    if (! registered_replica) {
      std::ostream *reg_os =
        cvm::proxy->output_stream(replicas_registry_file,
                                  std::ios::app);
      if (!reg_os) {
        return cvm::get_error();
      }
      *reg_os << replica_id << " " << replica_list_file << "\n";
      cvm::proxy->close_output_stream(replicas_registry_file);
    }
  }

  if (b_hills_traj) {
    if (!hills_traj_os) {
      hills_traj_os = cvm::proxy->output_stream(hills_traj_file_name());
      if (!hills_traj_os) return cvm::get_error();
    }
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
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
  if (this->comm != single_replica)
    os << "replicaID " << this->replica_id << "\n";
  return (colvarbias::get_state_params() + os.str());
}


std::ostream & colvarbias_meta::write_state_data(std::ostream& os)
{
  if (use_grids) {

    // this is a very good time to project hills, if you haven't done
    // it already!
    project_hills(new_hills_begin, hills.end(),
                  hills_energy,    hills_energy_gradients);
    new_hills_begin = hills.end();

    // write down the grids to the restart file
    os << "  hills_energy\n";
    hills_energy->write_restart(os);
    os << "  hills_energy_gradients\n";
    hills_energy_gradients->write_restart(os);
  }

  if ( (!use_grids) || keep_hills ) {
    // write all hills currently in memory
    for (std::list<hill>::const_iterator h = this->hills.begin();
         h != this->hills.end();
         h++) {
      os << *h;
    }
  } else {
    // write just those that are near the grid boundaries
    for (std::list<hill>::const_iterator h = this->hills_off_grid.begin();
         h != this->hills_off_grid.end();
         h++) {
      os << *h;
    }
  }
  if (ebmeta && ebmetaerror) {
    os << "ebmeta {\n";
    os << "    ebmetaNconst " << ebmeta_nconst << "\n"; 
    os << "    targetDist ";
    for (size_t j = 0; j < eff_error_points; j++) {
       os << target_dist_eff[j];
       os << " ";
    }
    os << "\n";
    os << "    gammaVec ";
    for (size_t j = 0; j < eff_error_points; j++) {
       os << gamma_vec[j];
       os << " ";
    }
    os << "\n";
    if (update_targets) { 
      os << "    targetError ";
      for (size_t j = 0; j < eff_error_points; j++) {
         size_t i = which_error_point[j];
         os << target_error->array_value(i);
         os << " ";
      }
      os << "\n";
      os << "    targetProb ";
      for (size_t j = 0; j < eff_error_points; j++) {
         os << target_prob[j];
         os << " ";
      }
      os << "\n"; 
    }
    os << "\n";
    os << "}\n";
  }

  return os;
}


int colvarbias_meta::write_state_to_replicas()
{
  if (comm != single_replica) {
    write_replica_state_file();
    // schedule to reread the state files of the other replicas (they
    // have also rewritten them)
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      (replicas[ir])->replica_state_file_in_sync = false;
    }
  }
  return COLVARS_OK;
}


int colvarbias_meta::write_output_files()
{
  if (dump_fes) {
    write_pmf();
  }
  return COLVARS_OK;
}


void colvarbias_meta::write_pmf()
{
  // allocate a new grid to store the pmf
  colvar_grid_scalar *pmf = new colvar_grid_scalar(*hills_energy);
  pmf->setup();

  if ((comm == single_replica) || (dump_replica_fes)) {
    // output the PMF from this instance or replica
    pmf->reset();
    pmf->add_grid(*hills_energy);
    cvm::real const max = pmf->maximum_value();
    pmf->add_constant(-1.0 * max);
    pmf->multiply_constant(-1.0);
    if (well_tempered) {
      cvm::real const well_temper_scale = (bias_temperature + cvm::temperature()) / bias_temperature;
      pmf->multiply_constant(well_temper_scale);
    }
    {
      std::string const fes_file_name(this->output_prefix +
                                      ((comm != single_replica) ? ".partial" : "") +
                                      (dump_fes_save ?
                                       "."+cvm::to_str(cvm::step_absolute()) : "") +
                                      ".pmf");
      cvm::proxy->backup_file(fes_file_name);
      std::ostream *fes_os = cvm::proxy->output_stream(fes_file_name);
      pmf->write_multicol(*fes_os);
      cvm::proxy->close_output_stream(fes_file_name);
    }
  }

  if (comm != single_replica) {
    // output the combined PMF from all replicas
    pmf->reset();
    pmf->add_grid(*hills_energy);
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      pmf->add_grid(*(replicas[ir]->hills_energy));
    }
    cvm::real const max = pmf->maximum_value();
    pmf->add_constant(-1.0 * max);
    pmf->multiply_constant(-1.0);
    if (well_tempered) {
      cvm::real const well_temper_scale = (bias_temperature + cvm::temperature()) / bias_temperature;
      pmf->multiply_constant(well_temper_scale);
    }
    std::string const fes_file_name(this->output_prefix +
                                    (dump_fes_save ?
                                     "."+cvm::to_str(cvm::step_absolute()) : "") +
                                    ".pmf");
    cvm::proxy->backup_file(fes_file_name);
    std::ostream *fes_os = cvm::proxy->output_stream(fes_file_name);
    pmf->write_multicol(*fes_os);
    cvm::proxy->close_output_stream(fes_file_name);
  }

  delete pmf;
}



int colvarbias_meta::write_replica_state_file()
{
  if (cvm::debug()) {
    cvm::log("Writing replica state file for bias \""+name+"\"\n");
  }
  // write down also the restart for the other replicas
  cvm::backup_file(replica_state_file.c_str());
  std::ostream *rep_state_os = cvm::proxy->output_stream(replica_state_file);
  if (rep_state_os == NULL) {
    cvm::error("Error: in opening file \""+
               replica_state_file+"\" for writing.\n", FILE_ERROR);
    return FILE_ERROR;
  }

  rep_state_os->setf(std::ios::scientific, std::ios::floatfield);

  if (!write_state(*rep_state_os)) {
    cvm::error("Error: in writing to file \""+
               replica_state_file+"\".\n", FILE_ERROR);
    cvm::proxy->close_output_stream(replica_state_file);
    return FILE_ERROR;
  }

  cvm::proxy->close_output_stream(replica_state_file);

  // rep_state_os.setf(std::ios::scientific, std::ios::floatfield);
  // rep_state_os << "\n"
  //              << "metadynamics {\n"
  //              << "  configuration {\n"
  //              << "    name " << this->name << "\n"
  //              << "    step " << cvm::step_absolute() << "\n";
  // if (this->comm != single_replica) {
  //   rep_state_os << "    replicaID " << this->replica_id << "\n";
  // }
  // rep_state_os << "  }\n\n";
  // rep_state_os << "  hills_energy\n";
  // rep_state_os << std::setprecision(cvm::cv_prec)
  //              << std::setw(cvm::cv_width);
  // hills_energy->write_restart(rep_state_os);
  // rep_state_os << "  hills_energy_gradients\n";
  // rep_state_os << std::setprecision(cvm::cv_prec)
  //              << std::setw(cvm::cv_width);
  // hills_energy_gradients->write_restart(rep_state_os);

  // if ( (!use_grids) || keep_hills ) {
  //   // write all hills currently in memory
  //   for (std::list<hill>::const_iterator h = this->hills.begin();
  //        h != this->hills.end();
  //        h++) {
  //     rep_state_os << *h;
  //   }
  // } else {
  //   // write just those that are near the grid boundaries
  //   for (std::list<hill>::const_iterator h = this->hills_off_grid.begin();
  //        h != this->hills_off_grid.end();
  //        h++) {
  //     rep_state_os << *h;
  //   }
  // }
  // rep_state_os << "}\n\n";
  // rep_state_os.close();

  // reopen the hills file
  cvm::proxy->close_output_stream(replica_hills_file);
  cvm::proxy->backup_file(replica_hills_file);
  replica_hills_os = cvm::proxy->output_stream(replica_hills_file);
  if (!replica_hills_os) return cvm::get_error();
  replica_hills_os->setf(std::ios::scientific, std::ios::floatfield);

  return COLVARS_OK;
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
  for (i = 0; i < widths.size(); i++) {
    os << " ";
    os << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width) << widths[i];
  }

  os << "  ";
  os << std::setprecision(cvm::en_prec)
     << std::setw(cvm::en_width) << W << "\n";

  return os.str();
}


std::ostream & operator << (std::ostream &os, colvarbias_meta::hill const &h)
{
  os.setf(std::ios::scientific, std::ios::floatfield);

  os << "hill {\n";
  os << "  step " << std::setw(cvm::it_width) << h.it << "\n";
  os << "  weight   "
     << std::setprecision(cvm::en_prec)
     << std::setw(cvm::en_width)
     << h.W << "\n";

  if (h.replica.size())
    os << "  replicaID  " << h.replica << "\n";

  size_t i;
  os << "  centers ";
  for (i = 0; i < (h.centers).size(); i++) {
    os << " "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << h.centers[i];
  }
  os << "\n";

  os << "  widths  ";
  for (i = 0; i < (h.widths).size(); i++) {
    os << " "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << h.widths[i];
  }
  os << "\n";

  os << "}\n";

  return os;
}


