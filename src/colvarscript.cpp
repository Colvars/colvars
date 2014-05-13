#include "colvarscript.h"

/// Parse config from file
int colvarscript::configfile (std::string const &filename) {
}

/// Parse config from string
int colvarscript::configstring (std::string const &config) {
}

/// delete object (type colvar or bias or any type)
/// A colvar may not be deleted if a bias depends on it (the bias should be deleted first)
int colvarscript::delete_obj (std::string const &object, std::string const &name) {
}

/// In analysis mode, go to specified frame (-1 for last frame)
/// returns new frame number, or COLVARSCRIPT_ERROR
int colvarscript::frame (int frame) {
}

/// In analysis mode, go to next frame, if it exists (-1 for last frame)
/// returns new frame number, or COLVARSCRIPT_ERROR
int colvarscript::nextframe () {
}

/// Recalculate all colvars and biases
int colvarscript::update () {
}

/// Delete every child object
int colvarscript::reset () {
}

/// Get colvar value
int colvarscript::get_value (std::string const &name, colvarvalue &value_out) {
  colvar *cv = cvm::colvar_by_name (name);
  if (cv == NULL) {
    error_message = "Colvar not found: " + name;
    return COLVARSCRIPT_ERROR;
  }
  value_out = cv->value();
  return COLVARSCRIPT_OK;
}

/// Get bias energy
int colvarscript::get_energy (std::string const &name, cvm::real &energy_out) {
  colvarbias *b = cvm::bias_by_name (name);
  if (b == NULL) {
    error_message = "Bias not found: " + name;
    return COLVARSCRIPT_ERROR;
  }
  energy_out = b->get_energy();
  return COLVARSCRIPT_OK;
}
