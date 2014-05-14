// -*- c++ -*-

#ifndef COLVARSCRIPT_H
#define COLVARSCRIPT_H

#include <string>
#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarbias.h"

#define COLVARSCRIPT_ERROR -1
#define COLVARSCRIPT_OK 0

class colvarscript  {

private:

public:

  inline colvarscript() {}
  inline ~colvarscript() {}

  /// If an error is return by one of the methods, it should set this to the error message
  std::string error_message;

  /// Parse the positional arguments of a script command
  int args (int argc, const char *argv[]);

  /// Parse config from file
  int configfile (std::string const &filename);

  /// Parse config from string
  int configstring (std::string const &config);

  /// delete object (type colvar or bias or any type)
  /// A colvar may not be deleted if a bias depends on it (the bias should be deleted first)
  int delete_obj (std::string const &object, std::string const &name);

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
};

#endif
