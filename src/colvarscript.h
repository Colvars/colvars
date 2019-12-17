// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARSCRIPT_H
//#define COLVARSCRIPT_H // Delay definition until later

#include <string>
#include <vector>
#include <map>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias.h"
#include "colvarproxy.h"


// Only these error values are part of the scripting interface
#define COLVARSCRIPT_ERROR -1
#define COLVARSCRIPT_OK 0


class colvarscript  {

private:

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// temporary until relevant functions become colvarscript methods
public:

  colvarproxy *proxy;
  colvarmodule *colvars;

  inline colvarscript() {} // no-argument construction forbidden

public:

  friend class colvarproxy;

  colvarscript(colvarproxy * p);
  inline ~colvarscript() {}

  /// If an error is caught by the proxy through fatal_error(), this is set to
  /// COLVARSCRIPT_ERROR
  int proxy_error;

  /// If an error is returned by one of the methods, it should set this to the
  /// error message
  std::string result;

  /// Run script command with given positional arguments (objects)
  int run(int objc, unsigned char *const objv[]);

  /// Set the return value of the script command to the given string
  inline void set_result_str(std::string const &s)
  {
    result = s;
  }

  /// Set the error message of the script interface to the given string
  inline void set_error_msg(std::string const &s)
  {
    result = s;
  }

  /// Build and return a short help
  std::string help_string(void) const;

  /// Use scripting language to get the string representation of an object
  inline char const *obj_to_str(unsigned char *const obj)
  {
    return cvm::proxy->script_obj_to_str(obj);
  }

  enum command {
    cv_help,
    cv_version,
    cv_config,
    cv_configfile,
    cv_getconfig,
    cv_reset,
    cv_resetindexgroups,
    cv_delete,
    cv_list,
    cv_list_biases,
    cv_load,
    cv_save,
    cv_update,
    cv_addenergy,
    cv_getenergy,
    cv_printframe,
    cv_printframelabels,
    cv_frame,
    cv_units,
    cv_colvar,
    // Should the following subcommands be listed in the enum if they are handled by the cv_colvar function?
    // do we want a mechanism to keep help strings for subcommands, or one big help string for the
    // top-level commands colvar and bias?
    cv_colvar_value,
    cv_colvar_update,
    cv_colvar_type,
    cv_colvar_delete,
    cv_colvar_addforce,
    cv_colvar_getappliedforce,
    cv_colvar_gettotalforce,
    cv_colvar_cvcflags,
    cv_colvar_getconfig,
    cv_colvar_get,
    cv_colvar_set,
    cv_bias,
    // Should the following subcommands be listed in the enum if they are handled by the cv_bias function?
    cv_bias_energy,
    cv_bias_update,
    cv_bias_delete,
    cv_bias_getconfig,
    cv_bias_get,
    cv_bias_set,
    cv_n_commands
  };

  /// Execute a script command
  inline int exec_command(command c,
                          void *pobj,
                          int objc, unsigned char * const *objv)
  {
    return (*(comm_fns[c]))(pobj, objc, objv);
  }

  /// Get help for a command (TODO reformat for each language?)
  inline std::string command_help(colvarscript::command c) const
  {
    return comm_help[c];
  }

  /// Clear all object results
  inline void clear_results()
  {
    result.clear();
  }

private:
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// temporary until relevant functions become colvarscript methods
public:

  /// Run subcommands on colvar
  int proc_colvar(colvar *cv, int argc, unsigned char *const argv[]);

  /// Run subcommands on bias
  int proc_bias(colvarbias *b, int argc, unsigned char *const argv[]);

  /// Run subcommands on base colvardeps object (colvar, bias, ...)
  int proc_features(colvardeps *obj,
                    int argc, unsigned char *const argv[]);

  /// Internal identifiers of command strings
  std::map<std::string, command> comm_str_map;

  /// Help strings for each command
  std::vector<std::string> comm_help;

  /// Number of arguments for each command
  std::vector<size_t> comm_n_args;

  /// Arguments for each command
  std::vector< std::vector<std::string> > comm_args;

  /// Implementations of each command
  std::vector<int (*)(void *, int, unsigned char * const *)> comm_fns;

};


/// Get a pointer to the main colvarscript object
inline static colvarscript *colvarscript_obj()
{
  return cvm::main()->proxy->script;
}

/// Get a pointer to the colvar object pointed to by pobj
inline static colvar *colvar_obj(void *pobj)
{
  return reinterpret_cast<colvar *>(pobj);
}

/// Get a pointer to the colvarbias object pointed to by pobj
inline static colvarbias *colvarbias_obj(void *pobj)
{
  return reinterpret_cast<colvarbias *>(pobj);
}


#define CVSCRIPT_COMM_FNAME(COMM) cvscript_ ## COMM

#define CVSCRIPT_COMM_PROTO(COMM)                                       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *, int, unsigned char *const *);

#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_PROTO(COMM)

#undef COLVARSCRIPT_H
#endif // #ifndef COLVARSCRIPT_H


#ifdef COLVARSCRIPT_CPP
#define CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)       \
  extern "C" int CVSCRIPT_COMM_FNAME(COMM)(void *pobj,                  \
                                           int objc,                    \
                                           unsigned char *const objv[]) \
  {                                                                     \
    colvarscript *script = colvarscript_obj();                          \
    script->clear_results();                                            \
    if (objc < 2+N_ARGS_MIN) /* "cv" and "COMM" are 1st and 2nd */ {    \
      script->set_error_msg("Missing arguments\n" +                    \
                             script->command_help(colvarscript::COMM)); \
      return COLVARSCRIPT_ERROR;                                        \
    }                                                                   \
    if (objc > 2+N_ARGS_MAX) {                                          \
      script->set_error_msg("Too many arguments\n" +                   \
                             script->command_help(colvarscript::COMM)); \
      return COLVARSCRIPT_ERROR;                                        \
    }                                                                   \
    FN_BODY;                                                            \
  }
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY) \
  CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)
#endif // #ifdef COLVARSCRIPT_CPP


#ifdef COLVARSCRIPT_INIT_FN
#define CVSCRIPT_COMM_INIT(COMM,HELP,ARGS) {                    \
    comm_str_map[#COMM] = COMM;                                 \
    comm_help[COMM] = HELP;                                     \
    comm_fns[COMM] = &(CVSCRIPT_COMM_FNAME(COMM));              \
  }
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_INIT(COMM,HELP,ARGS)
#endif


#if !defined(COLVARSCRIPT_H) || defined(COLVARSCRIPT_INIT_FN)
#define COLVARSCRIPT_H

#ifndef COLVARSCRIPT_INIT_FN
#ifdef __cplusplus
extern "C" {
#endif
#endif

  // Add optional arguments for command-specific help?
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
           { "conf (str) - Configuration string" },
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
           { "E (float) - Amount of energy to add" },
           cvm::main()->total_bias_energy +=
             strtod(script->obj_to_str(objv[2]), NULL);
           return COLVARS_OK;
           )

  CVSCRIPT(cv_getenergy,
           "Get the current Colvars energy",
           1, 1,
           { "E (float) - Store the energy in this variable" },
           double *energy = reinterpret_cast<double *>(objv[2]);
           *energy = cvm::main()->total_bias_energy;
           return COLVARS_OK;
           )

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


#ifndef COLVARSCRIPT_INIT_FN
#ifdef __cplusplus
} // extern "C"
#endif
#endif

#undef CVSCRIPT

#endif // #ifndef COLVARSCRIPT_H
