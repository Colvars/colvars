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


// Instantiate the body of all script commands

#define CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *pobj,                             \
                                int objc, unsigned char *const objv[])  \
  {                                                                     \
    colvarscript *script = colvarscript_obj();                          \
    script->clear_str_result();                                         \
    if (script->check_cmd_nargs<>(#COMM,                                \
                                  objc, N_ARGS_MIN, N_ARGS_MAX) !=      \
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
