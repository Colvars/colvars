// -*- c++ -*-

#include <vector>
#include <cstdlib>
#include <stdlib.h>
#include <string.h>

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
char const **cvscript_command_names()
{
  colvarscript *script = colvarscript_obj();
  return script->get_command_names();
}


extern "C"
char const *cvscript_help(char const *c)
{
  colvarscript *script = colvarscript_obj();
  return script->get_command_help(c).c_str();
}


// Get a pointer to the i-th argument (NULL if not given)
inline unsigned char *colvarscript_arg(int iarg,
                                       int objc, unsigned char *const objv[])
{
  // "cv" and "COMM" are 1st and 2nd
  return (2+iarg < objc) ? objv[2+iarg] : NULL;
}


// Check the argument count
inline int colvarscript_nargs_check(colvarscript *script,
                                    char const *comm_str,
                                    int objc,
                                    int n_args_min,
                                    int n_args_max)
{
  // "cv" and "COMM" are 1st and 2nd argument
  if (objc < 2+n_args_min) {
    script->set_result_str("Missing arguments\n" +
                           script->get_command_help(comm_str));
    return COLVARSCRIPT_ERROR;
  }
  if (objc > 2+n_args_max) {
    script->set_result_str("Too many arguments\n" +
                           script->get_command_help(comm_str));
    return COLVARSCRIPT_ERROR;
  }
  return COLVARSCRIPT_OK;
}


// Now instantiate the body of all script commands

#define CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *pobj,                             \
                                int objc, unsigned char *const objv[])  \
  {                                                                     \
    colvarscript *script = colvarscript_obj();                          \
    script->clear_str_result();                                         \
    if (colvarscript_nargs_check(script, #COMM,                         \
                                 objc, N_ARGS_MIN, N_ARGS_MAX) !=       \
        COLVARSCRIPT_OK) {                                              \
      return COLVARSCRIPT_ERROR;                                        \
    }                                                                   \
    FN_BODY;                                                            \
  }
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY) \
  CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)

#undef COLVARSCRIPT_COMMANDS_H
#include "colvarscript_commands.h"

#undef CVSCRIPT_COMM_FN
#undef CVSCRIPT
