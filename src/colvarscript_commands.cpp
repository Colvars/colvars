// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <vector>

#include "colvar.h"
#include "colvarbias.h"
#include "colvarproxy.h"
#include "colvardeps.h"
#include "colvarscript.h"
#include "colvarscript_commands.h"



extern "C"
int cvscript_n_commands()
{
  return static_cast<int>(colvarscript::cv_n_commands);
}


extern "C"
char const **cvscript_command_names(void *proxy_in)
{
  colvarproxy* proxy = (colvarproxy*)proxy_in;
  if (proxy) {
    if (proxy->script) {
      return proxy->script->get_command_names();
    }
  }
  return nullptr;
}


extern "C"
char const *cvscript_command_help(void *proxy_in, char const *c)
{
  colvarproxy* proxy = (colvarproxy*)proxy_in;
  if (proxy) {
    if (proxy->script) {
      return proxy->script->get_command_help(c);
    }
  }
  return nullptr;
}


extern "C"
char const *cvscript_command_rethelp(void *proxy_in, char const *c)
{
  colvarproxy* proxy = (colvarproxy*)proxy_in;
  if (proxy) {
    if (proxy->script) {
      return proxy->script->get_command_rethelp(c);
    }
  }
  return nullptr;
}


extern "C"
char const *cvscript_command_arghelp(void *proxy_in, char const *c, int i)
{
  colvarproxy* proxy = (colvarproxy*)proxy_in;
  if (proxy) {
    if (proxy->script) {
      return proxy->script->get_command_arghelp(c, i);
    }
  }
  return nullptr;
}


extern "C"
char const *cvscript_command_full_help(void *proxy_in, char const *c)
{
  colvarproxy* proxy = (colvarproxy*)proxy_in;
  if (proxy) {
    if (proxy->script) {
      return proxy->script->get_command_full_help(c);
    }
  }
  return nullptr;
}


extern "C"
int cvscript_command_n_args_min(void *proxy_in, char const *c)
{
  colvarproxy* proxy = (colvarproxy*)proxy_in;
  if (proxy) {
    if (proxy->script) {
      return proxy->script->get_command_n_args_min(c);
    }
  }
  return -1;
}


extern "C"
int cvscript_command_n_args_max(void *proxy_in, char const *c)
{
  colvarproxy* proxy = (colvarproxy*)proxy_in;
  if (proxy) {
    if (proxy->script) {
      return proxy->script->get_command_n_args_max(c);
    }
  }
  return -1;
}


// Instantiate the body of all script commands

#define CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *pobj,                             \
                                int objc, unsigned char *const objv[])  \
  {                                                                     \
    colvarscript *script = colvarscript_obj(pobj);                      \
    colvarmodule *cvmodule = script->module();                          \
    if (cvm::debug()) {                                            \
      cvmodule->log("Executing script function \""+std::string(#COMM)+"\""); \
    }                                                                   \
    script->clear_str_result();                                         \
    if (script->check_module_cmd_nargs(#COMM,                           \
                                       objc, N_ARGS_MIN, N_ARGS_MAX) != \
        COLVARSCRIPT_OK) {                                              \
      return COLVARSCRIPT_ERROR;                                        \
    }                                                                   \
    if (objc > 1) {                                                     \
      /* Silence unused parameter warning */                            \
      (void) objv[0];                                                   \
    }                                                                   \
    FN_BODY;                                                            \
  }
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY) \
  CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)

// Skips the colvar- and bias- specific commands
#define COLVARSCRIPT_COMMANDS_GLOBAL

#undef COLVARSCRIPT_COMMANDS_H
#include "colvarscript_commands.h"

#undef CVSCRIPT_COMM_FN
#undef CVSCRIPT
