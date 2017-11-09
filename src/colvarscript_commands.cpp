// -*- c++ -*-

#include <vector>
#include <cstdlib>
#include <stdlib.h>
#include <string.h>

#include "colvarproxy.h"
#include "colvardeps.h"
#include "colvarscript.h"
#include "colvarscript_commands.h"
#include "colvarscript_commands_util.h"



extern "C"
int run_colvarscript_command(int objc, unsigned char *const objv[])
{
  colvarscript *script = colvarscript_obj();
  if (!script) {
    cvm::error("Called run_colvarscript_command without a script object.\n",
               BUG_ERROR);
    return -1;
  }
  script->enter_interp_call(colvarscript::cv_text);
  int retval = script->run(objc, objv);
  script->exit_interp_call();
  return retval;
}


extern "C"
const char * get_colvarscript_result()
{
  colvarscript *script = colvarscript_obj();
  if (!script) {
    cvm::error("Called get_colvarscript_result without a script object.\n");
    return "";
  }
  return script->str_result().c_str();
}


extern "C"
int cvscript_n_commands()
{
  return static_cast<int>(colvarscript::cv_n_commands);
}


extern "C"
char const ** cvscript_command_names()
{
  colvarscript *script = colvarscript_obj();
  return script->get_command_names();
}


extern "C"
char const * cvscript_help(char const *c)
{
  colvarscript *script = colvarscript_obj();
  return script->get_command_help(c).c_str();
}


// extern "C"
// int (*cvscript_command(int c))(void *, int, unsigned char * const *)
// {
//   colvarscript *script = colvarscript_obj();
//   if ((c < 0) || (c > static_cast<int>(colvarscript::cv_n_commands))) {
//     cvm::error("Error: invalid script command index, "+cvm::to_str(c)+"\n",
//                BUG_ERROR);
//     return NULL;
//   }
//   colvarscript::command csc = static_cast<colvarscript::command>(c);
//   return script->get_command(csc);
// }


#define CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *pobj,                             \
                                int objc, unsigned char *const objv[])  \
  {                                                                     \
    colvarscript *script = colvarscript_obj();                          \
    script->clear_str_result();                                         \
    if (colvarscript_arg_check(script, #COMM,                           \
                               objc, N_ARGS_MIN, N_ARGS_MAX) !=         \
        COLVARSCRIPT_OK) {                                              \
      return COLVARSCRIPT_ERROR;                                        \
    }                                                                   \
    FN_BODY;                                                            \
  }
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY) \
  CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)

#undef COLVARSCRIPT_COMMANDS_H // disable include guard
#include "colvarscript_commands.h"

#undef CVSCRIPT_COMM_FN
#undef CVSCRIPT
