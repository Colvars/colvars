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

CVSCRIPT(cv_config,
         "Read configuration from the given string",
         1, 1,
         { "conf : str - Configuration string" },
         std::string const conf(script->obj_to_str(objv[2]));
         if (cvm::main()->read_config_string(conf) == COLVARS_OK) {
           return COLVARS_OK;
         }
         script->set_result_str("Error parsing configuration string");
         return COLVARSCRIPT_ERROR;
         )

CVSCRIPT(cv_getconfig,
         "Get the module's configuration string read so far",
         0, 0,
         { },
         script->set_result_str(cvm::main()->get_config());
         return COLVARS_OK;
         )

CVSCRIPT(cv_resetindexgroups,
         "Clear the index groups loaded so far, allowing to replace them",
         0, 0,
         { },
         cvm::main()->index_group_names.clear();
         cvm::main()->index_groups.clear();
         return COLVARS_OK;
         )

CVSCRIPT(cv_addenergy,
         "Add an energy to the MD engine",
         1, 1,
         { "E : float - Amount of energy to add" },
         cvm::main()->total_bias_energy +=
         strtod(script->obj_to_str(objv[2]), NULL);
         return COLVARSCRIPT_ERROR; // TODO Make this multi-language
         )

CVSCRIPT(cv_getenergy,
         "Get the current Colvars energy",
         1, 1,
         { "E : float - Copy the Colvars energy into this" },
         double const energy = cvm::main()->total_bias_energy;
         return script->set_result_real(energy,
                                        colvarscript_arg(1, objc, objv));
         )

CVSCRIPT(cv_getatomids,
         "Get the atom numbers of of the requested atoms",
         1, 1,
         { "a : int vector - Copy the atom IDs into this" },
         colvarproxy *proxy = cvm::main()->proxy;
         std::vector<int> const *a = proxy->modify_atom_ids();
         return script->set_result_int_vec(*a, colvarscript_arg(1, objc, objv));
         )

CVSCRIPT(cv_getatommasses,
         "Get the masses of the requested atoms",
         1, 1,
         { "m : float vector - Copy the masses into this" },
         colvarproxy *proxy = cvm::main()->proxy;
         std::vector<cvm::real> const *m = proxy->modify_atom_masses();
         return script->set_result_real_vec(*m,
                                            colvarscript_arg(1, objc, objv));
         )

CVSCRIPT(cv_getpositions,
         "Get the current positions of the requested atoms",
         1, 1,
         { "p : float vector - Copy the positions into this" },
         colvarproxy *proxy = cvm::main()->proxy;
         std::vector<cvm::rvector> const *p =
           proxy->modify_atom_positions();
         return script->set_result_rvector_vec(*p,
                                               colvarscript_arg(1, objc, objv));
         )

CVSCRIPT(cv_getappliedforces,
         "Get the current Colvars applied forces",
         1, 1,
         { "f : float vector - Copy the Colvars forces into this" },
         colvarproxy *proxy = cvm::main()->proxy;
         std::vector<cvm::rvector> const *f =
           proxy->modify_atom_new_colvar_forces();
         return script->set_result_rvector_vec(*f,
                                               colvarscript_arg(1, objc, objv));
         )

CVSCRIPT(cv_getatomgroupids,
         "Get the ID numbers of the requested atom groups",
         1, 1,
         { "a : int vector - Copy the atom group IDs into this" },
         colvarproxy *proxy = cvm::main()->proxy;
         std::vector<int> const *a = proxy->modify_atom_group_ids();
         return script->set_result_int_vec(*a, colvarscript_arg(1, objc, objv));
         )

CVSCRIPT(cv_getgrouppositions,
         "Get the current COM positions of the requested atom groups",
         1, 1,
         { "p : float vector - Copy the COM positions into this" },
         colvarproxy *proxy = cvm::main()->proxy;
         std::vector<cvm::rvector> const *p =
           proxy->modify_atom_group_positions();
         return script->set_result_rvector_vec(*p,
                                               colvarscript_arg(1, objc, objv));
         )

CVSCRIPT(cv_getgroupappliedforces,
         "Get the current COM forces applied to atom groups",
         1, 1,
         { "f : float vector - Copy the group forces into this" },
         colvarproxy *proxy = cvm::main()->proxy;
         std::vector<cvm::rvector> const *f =
           proxy->modify_atom_group_new_colvar_forces();
         return script->set_result_rvector_vec(*f,
                                               colvarscript_arg(1, objc, objv));
         )



#endif // #ifndef COLVARSCRIPT_COMMANDS_H
