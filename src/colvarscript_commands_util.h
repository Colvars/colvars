// -*- c++ -*-

#ifndef COLVARSCRIPT_COMMANDS_UTIL_H
#define COLVARSCRIPT_COMMANDS_UTIL_H

#include <vector>

#include "colvarproxy.h"
#include "colvarscript.h"
#include "colvarscript_commands.h"



// Get a pointer to the argument i-th (NULL if not given)
inline unsigned char *colvarscript_arg(int i,
                                       int objc, unsigned char *const objv[])
{
  // "cv" and "COMM" are 1st and 2nd
  return (2+i < objc) ? objv[2+i] : NULL;
}


/// Check the argument count
inline int colvarscript_arg_check(colvarscript *script,
                                  char const *comm_str,
                                  int objc,
                                  int n_args_min,
                                  int n_args_max)
{
  // "cv" and "COMM" are 1st and 2nd
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


inline std::vector<cvm::real>
pack_vector_rvectors(std::vector<cvm::rvector> const &v)
{
  std::vector<cvm::real> result(3*v.size(), 0.0);
  for (size_t i = 0; i < v.size(); i++) {
    result[3*i+0] = v[i].x;
    result[3*i+1] = v[i].y;
    result[3*i+2] = v[i].y;
  }
  return result;
}


#endif
