// -*- c++ -*-

#ifndef COLVARSCRIPT_H
#define COLVARSCRIPT_H

#include <string>
#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvarproxy.h"

#define COLVARSCRIPT_ERROR -1
#define COLVARSCRIPT_OK 0

class colvarscript  {

private:
  colvarproxy *proxy;
  colvarmodule *colvars;

  inline colvarscript() {} // no-argument construction forbidden

public:
 
  friend class colvarproxy;

  colvarscript(colvarproxy * p);
  inline ~colvarscript() {}

  /// If an error is caught by the proxy through fatal_error(), this is set to COLVARSCRIPT_ERROR
  int proxy_error;

  /// If an error is return by one of the methods, it should set this to the error message
  std::string result;

  /// Run script command with given positional arguments
  int run (int argc, char const *argv[]);

  /// Run subcommands on colvar
  int proc_colvar (int argc, char const *argv[]);

  /// Run subcommands on bias
  int proc_bias (int argc, char const *argv[]);

};

  
/*
  /// Parse config from file
  int configfile (std::string const &filename);

  /// Parse config from string
  int configstring (std::string const &config);

  /// delete object (type colvar or bias or any type)
  /// A colvar may not be deleted if a bias depends on it (the bias should be deleted first)
  int delete_obj (std::string const &object, std::string const &name);
  
  //Note: frame functions below are unneccessary in VMD as we rely on VMD's own current frame

  /// In analysis mode, go to specified frame (-1 for last frame)
  /// returns new frame number, or COLVARSCRIPT_ERROR
  int frame (int frame);

  /// In analysis mode, go to next frame, if it exists (-1 for last frame)
  /// returns new frame number, or COLVARSCRIPT_ERROR
  int nextframe ();

  /// Recalculate all colvars and biases
  int update ();

  /// Delete every child object
  int reset ();

  /// Get colvar value
  int get_value (std::string const &colvarname, colvarvalue &value_out);

  /// Get bias energy
  int get_energy (std::string const &biasname, cvm::real &energy_out);
*/

#endif
