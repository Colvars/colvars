// -*- c++ -*-

#ifndef COLVARSCRIPT_COMMANDS_H
#define COLVARSCRIPT_COMMANDS_H

// The following is a complete definition of the scripting API.

// The CVSCRIPT macro is used in four distinct contexts of use of this file:
// 1) Expand to the functions' prototypes (when included generically)
// 2) List colvarscript::command entries (when included in colvarscript.h)
// 3) Implement colvarscript::init() (when included in colvarscript.cpp)
// 4) Define the functions' bodies (when included in colvarscript_commands.cpp)


// Each command is created by an instance of the CVSCRIPT macro

// The arguments of the CVSCRIPT macro are:

// COMM = the id of the command (must be a member of colvarscript::command)

// HELP = a one-line description (C string literal) for the command

// N_ARGS_MIN = the lowest number of arguments allowed

// N_ARGS_MAX = the highest number of arguments allowed

// ARGS = an array of C string literals describing each parameter, with the
//        format "name : type[, optional] - description"

// FN_BODY = the implementation of the function; this should be a thin layer
//           over the existing classes; the "script" pointer to the
//           colvarscript object is already set by the CVSCRIPT_COMM_FN macro;
//           see also the functions in colvarscript_commands_util.h.

#ifndef CVSCRIPT_COMM_FNAME
#define CVSCRIPT_COMM_FNAME(COMM) cvscript_ ## COMM
#endif

// If CVSCRIPT is not defined, this file yields the function prototypes
#ifndef CVSCRIPT

#ifdef __cplusplus
#define CVSCRIPT_COMM_PROTO(COMM)                                       \
  extern "C" int CVSCRIPT_COMM_FNAME(COMM)(void *,                      \
                                           int, unsigned char *const *);
#else
#define CVSCRIPT_COMM_PROTO(COMM)                                       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *, int, unsigned char *const *);
#endif

#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_PROTO(COMM)

#endif


// TODO Add optional arguments for command-specific help?
CVSCRIPT(cv_help,
         "Print the help message",
         0, 0,
         {},
         script->set_result_str(script->help_string());
         return COLVARS_OK;
         )

  CVSCRIPT(cv_version,
           "Get the Colvars Module version number",
           0, 0,
           {},
           script->set_result_str(COLVARS_VERSION);
           return COLVARS_OK;
           )

CVSCRIPT(cv_config,
         "Read configuration from the given string",
         1, 1,
         { "conf : str - Configuration string" },
         std::string const conf(script->obj_to_str(objv[2]));
         if (cvm::main()->read_config_string(conf) == COLVARS_OK) {
           return COLVARS_OK;
         }
           script->set_error_msg("Error parsing configuration string");
           return COLVARSCRIPT_ERROR;
           )

  CVSCRIPT(cv_configfile,
           "Read configuration from a file",
           1, 1,
           {"conf_file (str) - Path to configuration file"},
            if (script->colvars->read_config_file(script->obj_to_str(objv[2])) == COLVARS_OK) {
              return COLVARS_OK;
            } else {
              script->set_error_msg("Error parsing configuration file");
              return COLVARSCRIPT_ERROR;
            }
           )

CVSCRIPT(cv_getconfig,
         "Get the module's configuration string read so far",
         0, 0,
         { },
         script->set_result_str(cvm::main()->get_config());
         return COLVARS_OK;
         )

  CVSCRIPT(cv_reset,
           "Delete all internal configuration",
           0, 0,
           {},
           script->colvars->reset();
           return COLVARS_OK;
           )

  CVSCRIPT(cv_resetindexgroups,
           "Clear the index groups loaded so far, allowing to replace them",
           0, 0,
           { },
           return cvm::main()->reset_index_groups();
           )

  CVSCRIPT(cv_delete,
           "Delete this Colvars module instance",
           0, 0,
           {},
            return script->proxy->request_deletion();
           )

  CVSCRIPT(cv_list,
           "Return a list of all variables",
           // For backward compatibility, accept argument "biases"
           0, 1,
           {},
            std::string res;
            if (objc == 2) {
              for (std::vector<colvar *>::iterator cvi = script->colvars->colvars.begin();
                  cvi != script->colvars->colvars.end();
                  ++cvi) {
                res += (cvi == script->colvars->colvars.begin() ? "" : " ") + (*cvi)->name;
              }
              script->set_result_str(res);
              return COLVARS_OK;
            } else if (!strcmp(script->obj_to_str(objv[2]), "biases")) {
              for (std::vector<colvarbias *>::iterator bi = script->colvars->biases.begin();
                  bi != script->colvars->biases.end();
                  ++bi) {
                res += (bi == script->colvars->biases.begin() ? "" : " ") + (*bi)->name;
              }
              script->set_result_str(res);
              return COLVARS_OK;
            } else {
              script->set_error_msg("Wrong arguments to command \"list\"\n" + script->help_string());
              return COLVARSCRIPT_ERROR;
            }
           )

  CVSCRIPT(cv_list_biases,
           "Return a list of all biases",
           0, 0,
           {},
            std::string res;
            for (std::vector<colvarbias *>::iterator bi = script->colvars->biases.begin();
                bi != script->colvars->biases.end();
                ++bi) {
              res += (bi == script->colvars->biases.begin() ? "" : " ") + (*bi)->name;
            }
            script->set_result_str(res);
            return COLVARS_OK;
           )

  CVSCRIPT(cv_load,
           "Load a state file (requires matching configuration)",
           1, 1,
           {"state_file (str) - Path to existing state file"},
            script->proxy->input_prefix() = script->obj_to_str(objv[2]);
            if (script->colvars->setup_input() == COLVARS_OK) {
              return COLVARS_OK;
            } else {
              script->set_error_msg("Error loading state file");
              return COLVARSCRIPT_ERROR;
            }
           )

  CVSCRIPT(cv_save,
           "Save state to a file",
           1, 1,
           {"state_file (str) - Path to state file"},
            script->proxy->output_prefix() = script->obj_to_str(objv[2]);
            int error = 0;
            error |= script->colvars->setup_output();
            error |= script->colvars->write_restart_file(script->colvars->output_prefix()+
                                                ".colvars.state");
            error |= script->colvars->write_output_files();
            return error ? COLVARSCRIPT_ERROR : COLVARS_OK;
           )

  CVSCRIPT(cv_update,
           "Recalculate colvars and biases",
           0, 0,
           {},
            int error_code = script->proxy->update_input();
            if (error_code) {
              script->set_error_msg("Error updating the Colvars module (input)");
              return error_code;
            }
            error_code |= script->colvars->calc();
            if (error_code) {
              script->set_error_msg("Error updating the Colvars module (calc)");
              return error_code;
            }
            error_code |= script->proxy->update_output();
            if (error_code) {
              script->set_error_msg("Error updating the Colvars module (output)");
            }
            return error_code;
           )

CVSCRIPT(cv_addenergy,
         "Add an energy to the MD engine",
         1, 1,
         { "E : float - Amount of energy to add" },
         cvm::main()->total_bias_energy +=
         strtod(script->obj_to_str(objv[2]), NULL);
           return COLVARS_OK;
//          return COLVARSCRIPT_ERROR; // TODO Make this multi-language
           )

  CVSCRIPT(cv_getenergy,
           "Get the current Colvars energy",
           1, 1,
           { "E (float) - Store the energy in this variable" },
           double *energy = reinterpret_cast<double *>(objv[2]);
           *energy = cvm::main()->total_bias_energy;
           return COLVARS_OK;
           )

// CVSCRIPT(cv_getenergy,
//          "Get the current Colvars energy",
//          1, 1,
//          { "E : float - Copy the Colvars energy into this" },
//          double const energy = cvm::main()->total_bias_energy;
//          return script->set_result_real(energy,
//                                         colvarscript_arg(1, objc, objv));
//          )

// CVSCRIPT(cv_getatomids,
//          "Get the atom numbers of of the requested atoms",
//          1, 1,
//          { "a : int vector - Copy the atom IDs into this" },
//          colvarproxy *proxy = cvm::main()->proxy;
//          std::vector<int> const *a = proxy->modify_atom_ids();
//          return script->set_result_int_vec(*a, colvarscript_arg(1, objc, objv));
//          )

// CVSCRIPT(cv_getatommasses,
//          "Get the masses of the requested atoms",
//          1, 1,
//          { "m : float vector - Copy the masses into this" },
//          colvarproxy *proxy = cvm::main()->proxy;
//          std::vector<cvm::real> const *m = proxy->modify_atom_masses();
//          return script->set_result_real_vec(*m,
//                                             colvarscript_arg(1, objc, objv));
//          )

// CVSCRIPT(cv_getpositions,
//          "Get the current positions of the requested atoms",
//          1, 1,
//          { "p : float vector - Copy the positions into this" },
//          colvarproxy *proxy = cvm::main()->proxy;
//          std::vector<cvm::rvector> const *p =
//            proxy->modify_atom_positions();
//          return script->set_result_rvector_vec(*p,
//                                                colvarscript_arg(1, objc, objv));
//          )

// CVSCRIPT(cv_getappliedforces,
//          "Get the current Colvars applied forces",
//          1, 1,
//          { "f : float vector - Copy the Colvars forces into this" },
//          colvarproxy *proxy = cvm::main()->proxy;
//          std::vector<cvm::rvector> const *f =
//            proxy->modify_atom_new_colvar_forces();
//          return script->set_result_rvector_vec(*f,
//                                                colvarscript_arg(1, objc, objv));
//          )

// CVSCRIPT(cv_getatomgroupids,
//          "Get the ID numbers of the requested atom groups",
//          1, 1,
//          { "a : int vector - Copy the atom group IDs into this" },
//          colvarproxy *proxy = cvm::main()->proxy;
//          std::vector<int> const *a = proxy->modify_atom_group_ids();
//          return script->set_result_int_vec(*a, colvarscript_arg(1, objc, objv));
//          )

// CVSCRIPT(cv_getgrouppositions,
//          "Get the current COM positions of the requested atom groups",
//          1, 1,
//          { "p : float vector - Copy the COM positions into this" },
//          colvarproxy *proxy = cvm::main()->proxy;
//          std::vector<cvm::rvector> const *p =
//            proxy->modify_atom_group_positions();
//          return script->set_result_rvector_vec(*p,
//                                                colvarscript_arg(1, objc, objv));
//          )

// CVSCRIPT(cv_getgroupappliedforces,
//          "Get the current COM forces applied to atom groups",
//          1, 1,
//          { "f : float vector - Copy the group forces into this" },
//          colvarproxy *proxy = cvm::main()->proxy;
//          std::vector<cvm::rvector> const *f =
//            proxy->modify_atom_group_new_colvar_forces();
//          return script->set_result_rvector_vec(*f,
//                                                colvarscript_arg(1, objc, objv));
//          )


  CVSCRIPT(cv_units,
           "Get the current Colvars unit system",
           0, 1,
           { },
           if (objc < 3) {
            script->set_result_str(cvm::proxy->units);
            return COLVARS_OK;
           } else {
            return cvm::proxy->set_unit_system(script->obj_to_str(objv[2]) , false);
           }
           )

  CVSCRIPT(cv_printframelabels,
           "Print the labels that would be written to colvars.traj",
           0, 0,
           { },
            std::ostringstream os;
            script->colvars->write_traj_label(os);
            script->set_result_str(os.str());
            return COLVARS_OK;
           )

  CVSCRIPT(cv_printframe,
           "Print the values that would be written to colvars.traj",
           0, 0,
           { },
            std::ostringstream os;
            script->colvars->write_traj(os);
            script->set_result_str(os.str());
            return COLVARS_OK;
           )

  CVSCRIPT(cv_frame,
           "Get or set current frame number",
           0, 1,
           { },
            if (objc == 2) {
              long int f;
              int error = script->proxy->get_frame(f);
              if (error == COLVARS_OK) {
                script->set_result_str(cvm::to_str(f));
                return COLVARS_OK;
              } else {
                script->set_error_msg("Frame number is not available");
                return COLVARSCRIPT_ERROR;
              }
            } else if (objc == 3) {
              // Failure of this function does not trigger an error, but
              // returns nonzero, to let scripts detect available frames
              int error = script->proxy->set_frame(strtol(script->obj_to_str(objv[2]), NULL, 10));
              script->set_result_str(cvm::to_str(error == COLVARS_OK ? 0 : -1));
              return COLVARS_OK;
            }
           )

  CVSCRIPT(cv_colvar,
           "Access colvar-specific commands: ",
           2, 4,
           {"name (str) - colvar name"},
            std::string const name(script->obj_to_str(objv[2]));
            colvar *cv = cvm::colvar_by_name(name);
            if (cv == NULL) {
              script->set_error_msg("Colvar not found: " + name);
              return COLVARSCRIPT_ERROR;
            }
            return script->proc_colvar(cv, objc-1, &(objv[1]));
           )

  CVSCRIPT(cv_bias,
           "Access bias-specific commands: ",
           2, 4,
           {"name (str) - bias name"},
            std::string const name(script->obj_to_str(objv[2]));
            colvarbias *cvb = cvm::bias_by_name(name);
            if (cvb == NULL) {
              script->set_error_msg("Bias not found: " + name);
              return COLVARSCRIPT_ERROR;
            }
            return script->proc_bias(cvb, objc-1, &(objv[1]));
           )




#endif // #ifndef COLVARSCRIPT_COMMANDS_H
