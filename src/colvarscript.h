// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARSCRIPT_H
#define COLVARSCRIPT_H

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

  colvarproxy *proxy;
  colvarmodule *colvars;

  inline colvarscript() {} // no-argument construction forbidden

public:

  friend class colvarproxy;

  colvarscript(colvarproxy *p);

  ~colvarscript();

  /// If an error is caught by the proxy through fatal_error(), this is set to
  /// COLVARSCRIPT_ERROR
  int proxy_error;

  /// If an error is returned by one of the methods, it should set this to the
  /// error message
  std::string result;

  /// Run script command with given positional arguments (objects)
  int run(int objc, unsigned char *const objv[]);

  /// Get the string result (if given)
  inline std::string const &str_result()
  {
    return result;
  }

  /// Set the return value to the given string
  int set_result_str(std::string const &s);

  /// Clear the string result
  int clear_str_result();

  /// Build and return a short help
  std::string help_string(void) const;

  /// Commands available
  enum command {
#define CVSCRIPT_ENUM_COMM(COMM) COMM,
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_ENUM_COMM(COMM)
#ifdef COLVARSCRIPT_COMMANDS_H
#undef COLVARSCRIPT_COMMANDS_H
#endif
#include "colvarscript_commands.h"
#undef COLVARSCRIPT_COMMANDS_H
#undef CVSCRIPT
#undef CVSCRIPT_ENUM_COMM
    cv_n_commands
  };


  /// Use scripting language to get the string representation of an object
  inline char const *obj_to_str(unsigned char *const obj)
  {
    return cvm::proxy->script_obj_to_str(obj);
  }
  
  /// Get names of all commands
  inline char const **get_command_names() const
  {
    return comm_names;
  }

  /// Get help string for a command
  std::string get_command_help(char const *c);

  /// Set error code for unsupported script operation
  int unsupported_op();

private:

  /// Set up all script API functions
  int init_commands();

  /// Set up a single script API function
  int init_command(colvarscript::command const &comm,
                   char const *name, char const *help,
                   int n_args_min, int n_args_max, char const **arghelp,
                   int (*fn)(void *, int, unsigned char * const *));

  /// Execute a script command
  inline int exec_command(command c,
                          void *pobj,
                          int objc, unsigned char * const *objv)
  {
    return (*(comm_fns[c]))(pobj, objc, objv);
  }

  /// Run subcommands on colvar
  int proc_colvar(colvar *cv, int argc, unsigned char *const argv[]);

  /// Run subcommands on bias
  int proc_bias(colvarbias *b, int argc, unsigned char *const argv[]);

  /// Run subcommands on base colvardeps object (colvar, bias, ...)
  int proc_features(colvardeps *obj,
                    int argc, unsigned char *const argv[]);

  /// Internal identifiers of command strings
  std::map<std::string, command> comm_str_map;

  /// Inverse of comm_str_map (to be exported outside this class)
  char const **comm_names;

  /// Help strings for each command
  std::vector<std::string> comm_help;

  /// Minimum number of arguments for each command
  std::vector<size_t> comm_n_args_min;

  /// Maximum number of arguments for each command
  std::vector<size_t> comm_n_args_max;

  /// Help strings for each command argument
  std::vector< std::vector<std::string> > comm_arghelp;

  /// Implementations of each command
  std::vector<int (*)(void *, int, unsigned char * const *)> comm_fns;

  /// Get a pointer to the implementation of the given command
  inline int (*get_comm_fn(std::string const &cmd_key))(void *,
                                                        int,
                                                        unsigned char * const *)
  {
    if (comm_str_map.count(cmd_key) > 0) {
      return comm_fns[comm_str_map[cmd_key]];
    }
    return NULL;
  }

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


extern "C" {

  /// Generic wrapper for string-based scripting
  int run_colvarscript_command(int objc, unsigned char *const objv[]);

  /// Get the string result of a script call
  const char * get_colvarscript_result();

}


#endif // #ifndef COLVARSCRIPT_H
