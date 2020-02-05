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

  colvarproxy *proxy_;
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

  /// Set the error message of the script interface to the given string
  inline void set_error_msg(std::string const &s)
  {
    result = s;
  }

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

  /// Type of object handling a script command
  enum Object_type {
    use_module,
    use_colvar,
    use_bias
  };

  /// Get a pointer to the i-th argument of the command (NULL if not given)
  template<Object_type T = use_module>
  unsigned char *get_cmd_arg(int iarg, int objc, unsigned char *const objv[]);

  /// Check the argument count of the command
  template<Object_type T = use_module>
  int check_cmd_nargs(char const *cmd, int objc,
                      int n_args_min, int n_args_max);

  /// Number of positional arguments to shift for each object type
  template<colvarscript::Object_type T>
  int cmd_arg_shift();

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

  /// Pointer to the Colvars main object
  inline colvarmodule *module()
  {
    return this->colvars;
  }

  /// Pointer to the colvarproxy object (interface with host engine)
  inline colvarproxy *proxy()
  {
    return this->proxy_;
  }

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



template<colvarscript::Object_type T>
int colvarscript::cmd_arg_shift()
{
  int shift = 0;
  if (T == use_module) {
    // "cv" and "COMMAND" are 1st and 2nd argument, and shift is equal to 2
    shift = 2;
  } else if (T == use_colvar) {
    // Same as above with additional arguments "colvar" and "NAME"
    shift = 4;
  } else if (T == use_bias) {
    shift = 4;
  }
  return shift;
}


template<colvarscript::Object_type T>
unsigned char *colvarscript::get_cmd_arg(int iarg,
                                         int objc,
                                         unsigned char *const objv[])
{
  int const shift = cmd_arg_shift<T>();
  return (shift+iarg < objc) ? objv[shift+iarg] : NULL;
}


template<colvarscript::Object_type T>
int colvarscript::check_cmd_nargs(char const *cmd,
                                  int objc,
                                  int n_args_min,
                                  int n_args_max)
{
  int const shift = cmd_arg_shift<T>();
  if (objc < shift+n_args_min) {
    set_result_str("Missing arguments\n" +
                   get_command_help(cmd));
    return COLVARSCRIPT_ERROR;
  }
  if (objc > shift+n_args_max) {
    set_result_str("Too many arguments\n" +
                   get_command_help(cmd));
    return COLVARSCRIPT_ERROR;
  }
  return COLVARSCRIPT_OK;
}


extern "C" {

  /// Generic wrapper for string-based scripting
  int run_colvarscript_command(int objc, unsigned char *const objv[]);

  /// Get the string result of a script call
  const char * get_colvarscript_result();

}


#endif // #ifndef COLVARSCRIPT_H
