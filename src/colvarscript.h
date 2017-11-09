// -*- c++ -*-

#ifndef COLVARSCRIPT_H
#define COLVARSCRIPT_H

#include <string>
#include <vector>
#include <map>

#if defined(NAMD_TCL) || defined(VMDTCL)
#define COLVARS_TCL
#include <tcl.h>
#endif

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

  /// Set up all script API functions
  int init_commands();

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

  /// Set the return value of the script command to the given string
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
  inline char const *obj_to_str(unsigned char *obj)
  {
    return cvm::proxy->script_obj_to_str(obj);
  }

  /// Execute a script command (TODO: deprecate after conversion)
  inline int exec_command(command c,
                          void *pobj,
                          int objc, unsigned char * const *objv)
  {
    return (*(comm_fns[c]))(pobj, objc, objv);
  }

  /// Get names of all a command
  inline char const ** get_command_names()
  {
    return comm_names;
  }

  /// Get help string for a command
  std::string get_command_help(char const *c);

  // /// Get a script command's function pointer
  // int (*get_command(char const *c))(void *, int,
  //                                   unsigned char * const *);

  /// Track calls from interpreters of different types
  enum interp_type {
    cv_text, // Use string representations (fallback for all other cases)
    cv_tcl, // Use Tcl object API (ref counting included)
    cv_nointerp // Use pointers to C variables/arrays (no ref count)
  };

  /// Enter a call to the interpreter-agnostic code
  inline void enter_interp_call(interp_type t)
  {
    clear_str_result();
    interp_type_list.push_back(t);
    interp_obj_list.push_back(NULL);
    interp_obj_size_list.push_back(0);
  }

  /// Exit a call to the interpreter-agnostic code
  inline void exit_interp_call()
  {
    interp_type_list.pop_back();
    interp_obj_list.pop_back();
    interp_obj_size_list.pop_back();
  }

  /// Pointer to the object returned by the current script call
  inline unsigned char *obj_result()
  {
    return interp_obj_list.back();
  }

  /// Size (in bytes) of obj_result(); -1 if using an external library
  inline int obj_result_size()
  {
    return interp_obj_size_list.back();
  }

  /// Modify obj_result()
  inline unsigned char **modify_obj_result()
  {
    return &(interp_obj_list.back());
  }

  /// Modify obj_result_size()
  inline int *modify_obj_result_size()
  {
    return &(interp_obj_size_list.back());
  }

  /// Copy the number x into obj if not NULL, or into obj_result() otherwise
  int set_result_real(cvm::real x, unsigned char *obj);

  /// Copy the vector x into obj if not NULL, or into obj_result() otherwise
  int set_result_real_vec(std::vector<cvm::real> const &x, unsigned char *obj);

  /// Copy the vector x into obj if not NULL, or into obj_result() otherwise
  int set_result_int_vec(std::vector<int> const &x, unsigned char *obj);

  /// Set error for unsupported script operation
  int unsupported_op();

private:

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

  /// Track calls from interpreters of different types
  std::list<interp_type> interp_type_list;

  /// Pointer to the result of each call
  std::list<unsigned char *> interp_obj_list;

  /// Size of the object return by each call
  std::list<int> interp_obj_size_list;
};


// Get a pointer to the main colvarscript object
inline colvarscript *colvarscript_obj()
{
  colvarmodule *cv = cvm::main();
  if (cv) {
    return cvm::main()->proxy->script;
  }
  return NULL;
}


// TODO Move this to colvarscript_capi.h
extern "C" {

#if defined(COLVARS_TCL)

  /// Tcl wrapper
  int tcl_run_colvarscript_command(ClientData clientData,
                                   Tcl_Interp *interp_in,
                                   int objc, Tcl_Obj *const objv[]);

#endif // #if defined(COLVARS_TCL)

  /// Generic wrapper for string-based scripting
  int run_colvarscript_command(int objc, unsigned char *const objv[]);

  /// Get the string result of a script call
  const char * get_colvarscript_result();

  /// Get the number of colvarscript commands
  int cvscript_n_commands();

  /// Get the names of all commands (array of strings)
  char const ** cvscript_commands();

  /// Get the help string of the given command (split using newlines)
  char const * cvscript_help(char const *cmd);

  // /// Get the function pointer of the given command
  // int (*cvscript_command(int c))(void *, int, unsigned char * const *);

}


#endif // #ifndef COLVARSCRIPT_H
