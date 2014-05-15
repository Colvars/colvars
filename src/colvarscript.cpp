#include "colvarscript.h"

/// Run method based on given arguments
int colvarscript::run (int argc, char *argv[]) {

  result = "";

  if (cvm::debug()) {
    cvm::log ("Called script run with " + cvm::to_str(argc) + " args");
    for (int i = 0; i < argc; i++) { cvm::log (argv[i]); }
  }

  if (argc == 1) {
    // TODO display usage statement
    // need console output (see proxy)
    // or put it in error message string
    return COLVARSCRIPT_OK;
  }

  std::string cmd = argv[1];

  if (argc == 2) {
    if (cmd == "reset") {
      /// Delete every child object
      // Implementation postponed until delayed initialization is in place
      return COLVARSCRIPT_OK;
    }
    
    if (cmd == "update") { 
        /// Recalculate all colvars and biases
      return COLVARSCRIPT_OK;
    }

    result = "Syntax error";
    return COLVARSCRIPT_ERROR;
  }

  if (argc == 3) {

    /// Parse config from file
    //int colvarscript::configfile (std::string const &filename) {
      // Implementation postponed until delayed initialization is in place

    /// Parse config from string
    if (cmd == "config") {
      std::string conf = argv[2];
      // Partial implementation: we are not reading global options
      colvars->init_colvars (conf);
      colvars->init_biases (conf);
      return COLVARSCRIPT_OK;
    }

    result = "Syntax error";
    return COLVARSCRIPT_ERROR;
  }

  if (cmd == "colvar") {
    std::string name = argv[2];
    colvar *cv = cvm::colvar_by_name (name);
    if (cv == NULL) {
      result = "Colvar not found: " + name;
      return COLVARSCRIPT_ERROR;
    }

    if (argc == 4) {
      result = "Syntax error";
      return COLVARSCRIPT_ERROR;
    }

    if (argc >= 5) {
      std::string subcmd = argv[3];
      std::string param = argv[4];

      if (subcmd == "delete") {
        if (cv->biases.size() > 0) {
          result = "Cannot delete a colvar currently used by biases, delete those biases first";
          return COLVARSCRIPT_ERROR;
        }
        //cvm::delete_colvar (cv); // TODO, implement, print new legend line in traj file
        result = "Function is not impemented yet";
        return COLVARSCRIPT_ERROR;
      }

      result = "Syntax error";
      return COLVARSCRIPT_ERROR;
    }
  }
  if (cmd == "bias") {
    std::string name = argv[2];
    colvarbias *b = cvm::bias_by_name (name);
    if (b == NULL) {
      result = "Bias not found: " + name;
      return COLVARSCRIPT_ERROR;
    }

    if (argc == 4) {
      result = "Syntax error";
      return COLVARSCRIPT_ERROR;
    }

    if (argc >= 5) {
      std::string subcmd = argv[3];
      std::string param = argv[4];

      if (subcmd == "delete") {
        //cvm::delete_bias (b); // TODO implement
        result = "Function is not impemented yet";
        return COLVARSCRIPT_ERROR;
      }
      result = "Syntax error";
      return COLVARSCRIPT_ERROR;
    }
  }
  result = "Syntax error";
  return COLVARSCRIPT_ERROR;
}

/*
/// Get colvar value
int colvarscript::get_value (std::string const &name, colvarvalue &value_out) {
  colvar *cv = cvm::colvar_by_name (name);
  if (cv == NULL) {
    result = "Colvar not found: " + name;
    return COLVARSCRIPT_ERROR;
  }
  value_out = cv->value();
  return COLVARSCRIPT_OK;
}

/// Get bias energy
int colvarscript::get_energy (std::string const &name, cvm::real &energy_out) {
  colvarbias *b = cvm::bias_by_name (name);
  if (b == NULL) {
    result = "Bias not found: " + name;
    return COLVARSCRIPT_ERROR;
  }
  energy_out = b->get_energy();
  return COLVARSCRIPT_OK;
}
*/
