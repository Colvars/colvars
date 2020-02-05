// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


CVSCRIPT(bias_bin,
         "Get the current grid bin index (1D ABF only for now)",
         0, 0,
         "",
         script->set_result_str(cvm::to_str(this_bias->current_bin()));
         return COLVARS_OK;
         )

CVSCRIPT(bias_bincount,
         "Get the number of samples at the given grid bin (1D ABF only for now)",
         0, 1,
         "index : integer - Grid index; defaults to current bin",
         int index = this_bias->current_bin();
         char const *indexarg =
           script->obj_to_str(script->get_cmd_arg<colvarscript::use_bias>(0, objc, objv));
         if (indexarg) {
           std::string const param(indexarg);
           if (!(std::istringstream(param) >> index)) {
             script->set_error_msg("bincount: error parsing bin index");
             return COLVARSCRIPT_ERROR;
           }
         }
         script->set_result_str(cvm::to_str(this_bias->bin_count(index)));
         return COLVARS_OK;
         )

CVSCRIPT(bias_binnum,
         "Get the total number of grid points of this bias (1D ABF only for now)",
         0, 0,
         "",
         int r = this_bias->bin_num();
         if (r < 0) {
           script->set_error_msg("Error: calling bin_num() for bias " +
                                 this_bias->name);
           return COLVARSCRIPT_ERROR;
         }
         script->set_result_str(cvm::to_str(r));
         return COLVARS_OK;
         )

CVSCRIPT(bias_delete,
         "Delete this bias",
         0, 0,
         "",
         delete this_bias;
         return COLVARS_OK;
         )

CVSCRIPT(bias_energy,
         "Get the current energy of this bias",
         0, 0,
         "",
         script->set_result_str(cvm::to_str(this_bias->get_energy()));
         return COLVARS_OK;
         )

CVSCRIPT(bias_get,
         "Get the value of the given feature for this bias",
         1, 1,
         "feature : string - Name of the feature",
         return script->proc_features(this_bias, objc, objv);
         )

CVSCRIPT(bias_getconfig,
         "Return the configuration string of this bias",
         0, 0,
         "",
         script->set_result_str(this_bias->get_config());
         return COLVARS_OK;
         )

CVSCRIPT(bias_set,
         "Set the given feature of this bias to a new value",
         2, 2,
         "feature : string - Name of the feature\n"
         "value : string - String representation of the new feature value",
         return script->proc_features(this_bias, objc, objv);
         )

CVSCRIPT(bias_share,
         "Share bias information with other replicas (multiple-walker scheme)",
         0, 0,
         "",
         if (this_bias->replica_share() != COLVARS_OK) {
           script->set_error_msg("Error: calling replica_share() for bias " +
                                 this_bias->name);
           return COLVARSCRIPT_ERROR;
         }
         return COLVARS_OK;
         )

CVSCRIPT(bias_state,
         "Return a string representation of the feature state of this bias",
         0, 0,
         "",
         return script->proc_features(this_bias, objc, objv);
         )

CVSCRIPT(bias_update,
         "Recompute this bias and return its up-to-date energy",
         0, 0,
         "",
         this_bias->update();
         script->set_result_str(cvm::to_str(this_bias->get_energy()));
         return COLVARS_OK;
         )
