// -*- c++ -*-

#include <cstdlib>
#include "colvarscript.h"


colvarscript::colvarscript (colvarproxy *p)
 : proxy (p)
{
  colvars = proxy->colvars;
}

/// Run method based on given arguments
int colvarscript::run (int argc, char const *argv[]) {

  result = "";

  if (cvm::debug()) {
    cvm::log ("Called script run with " + cvm::to_str(argc) + " args");
    for (int i = 0; i < argc; i++) { cvm::log (argv[i]); }
  }

  if (argc < 2) {
    result = "usage: "+std::string (argv[0])+"<subcommand> [args...]\n\
\n\
Initialize or deinitialize the module:\n\
  configfile <file name>      -- read configuration from a file\n\
  config <string>             -- read configuration from the given string\n\
  reset                       -- delete all internal configuration\n\
\n\
Input and output:\n\
  list                        -- return a list of all variables\n\
  list biases                 -- return a list of all biases\n\
  load <file name>            -- load a state file (requires configuration)\n\
  update                      -- recalculate colvars and biases based\n\
  printframe                  -- return a summary of the current frame\n\
  printframelabels            -- return labels to annotate printframe's output\n";

  if (proxy->frame() != COLVARS_NOT_IMPLEMENTED) {
      result += "\
  frame                       -- return current frame number\n\
  frame <new_frame>           -- return current frame number\n";
  }

  result += "\n\
Access collective variables:\n\
  colvar <name> value         -- return the current value of the colvar <name>\n\
  colvar <name> update        -- recalculate the colvar <name>\n\
  colvar <name> delete        -- delete the colvar <name>\n\
\n\
Access biases and algorithms:\n\
  bias <name> energy          -- return the current energy of the bias <name>\n\
  bias <name> update          -- recalculate the bias <name>\n\
  bias <name> delete          -- delete the bias <name>\n\
\n\
";
    return COLVARSCRIPT_OK;
  }

  std::string cmd = argv[1];
  
  if (cmd == "colvar") {
    return proc_colvar (argc-1, &(argv[1])); 
  }

  if (cmd == "bias") {
    return proc_bias (argc-1, &(argv[1])); 
  }

  if (cmd == "reset") {
    /// Delete every child object
    colvars->reset();
    return COLVARSCRIPT_OK;
  }
  
  if (cmd == "update") { 
    colvars->calc();
    return COLVARSCRIPT_OK;
  }

  if (cmd == "list") {
    if (argc == 2) {
      for (std::vector<colvar *>::iterator cvi = colvars->colvars.begin();
           cvi != colvars->colvars.end();
           cvi++) {
        result += (cvi == colvars->colvars.begin() ? "" : " ") + (*cvi)->name;
      }
      return COLVARSCRIPT_OK;
    } else if (argc == 3 && !strcmp(argv[2], "biases")) {
      for (std::vector<colvarbias *>::iterator bi = colvars->biases.begin();
           bi != colvars->biases.end();
           bi++) {
        result += (bi == colvars->biases.begin() ? "" : " ") + (*bi)->name;
      }
      return COLVARSCRIPT_OK;
    } else {
      result = "Wrong arguments to command \"list\"";
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Parse config from file
  if (cmd == "configfile") {
    if (argc < 3) {
      result = "Missing arguments";
      return COLVARSCRIPT_ERROR;
    }
    colvars->config_file (argv[2]);
    return COLVARSCRIPT_OK;
  }

  /// Parse config from string
  if (cmd == "config") {
    if (argc < 3) {
      result = "Missing arguments";
      return COLVARSCRIPT_ERROR;
    }
    std::string conf = argv[2];
    colvars->config_string (conf);
    return COLVARSCRIPT_OK;
  }

  /// Load an input state file
  if (cmd == "load") {
    if (argc < 3) {
      result = "Missing arguments";
      return COLVARSCRIPT_ERROR;
    }
    proxy->input_prefix_str = argv[2];
    colvars->setup_input();
    return COLVARSCRIPT_OK;
  }

  /// TODO Write an output state file? (Useful for testing)

  /// Print the values that would go on colvars.traj
  if (cmd == "printframelabels") {
    std::ostringstream os;
    colvars->write_traj_label (os);
    result = os.str();
    return COLVARSCRIPT_OK;
  }
  if (cmd == "printframe") {
    std::ostringstream os;
    colvars->write_traj (os);
    result = os.str();
    return COLVARSCRIPT_OK;
  }

  if (cmd == "frame") {
    if (argc == 2) {
      int f = proxy->frame();
      if (f >= 0) {
        result = cvm::to_str (f);
        return COLVARSCRIPT_OK;
      } else {
        result = "Frame number is not available";
        return COLVARSCRIPT_ERROR;
      }
    } else if (argc == 3) {
      // Failure of this function does not trigger an error, but
      // returns the plain result to let scripts detect available frames
      long int f = proxy->frame(strtol(argv[2], NULL, 10));
      colvars->it = proxy->frame();
      result = cvm::to_str (f);
      return COLVARSCRIPT_OK;
    } else {
      result = "Wrong arguments to command \"frame\"";
      return COLVARSCRIPT_ERROR;
    }
  }

  result = "Syntax error";
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_colvar (int argc, char const *argv[]) {
  std::string name = argv[1];
  colvar *cv = cvm::colvar_by_name (name);
  if (cv == NULL) {
    result = "Colvar not found: " + name;
    return COLVARSCRIPT_ERROR;
  }
  if (argc < 3) {
    result = "Missing parameters";
    return COLVARSCRIPT_ERROR;
  }
  std::string subcmd = argv[2];

  if (subcmd == "value") {
    result = cvm::to_str(cv->value(), 0, 14);
    return COLVARSCRIPT_OK;
  }

  if (subcmd == "update") {
    // note: this is not sufficient for a newly created colvar
    // as atom coordinates may not be properly loaded
    // a full CVM update is required
    // or otherwise updating atom positions
    cv->update();
    result = cvm::to_str(cv->value(), 0, 14);
    return COLVARSCRIPT_OK;
  }

  if (subcmd == "delete") {
    if (cv->biases.size() > 0) {
      result = "Cannot delete a colvar currently used by biases, delete those biases first";
      return COLVARSCRIPT_ERROR;
    }
    // colvar destructor is tasked with the cleanup
    delete cv;
    // TODO this could be done by the destructors
    colvars->write_traj_label (colvars->cv_traj_os);
    return COLVARSCRIPT_OK;
  }

  result = "Syntax error";
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_bias (int argc, char const *argv[]) {
  std::string name = argv[1];
  colvarbias *b = cvm::bias_by_name (name);
  if (b == NULL) {
    result = "Bias not found: " + name;
    return COLVARSCRIPT_ERROR;
  }

  if (argc < 3) {
    result = "Missing parameters";
    return COLVARSCRIPT_ERROR;
  }
  std::string subcmd = argv[2];

  if (subcmd == "energy") {
    result = cvm::to_str(b->get_energy());
    return COLVARSCRIPT_OK;
  }

  if (subcmd == "update") {
    b->update();
    result = cvm::to_str(b->get_energy());
    return COLVARSCRIPT_OK;
  }

  if (subcmd == "delete") {
    // the bias destructor takes care of the cleanup at cvm level
    delete b;
    // TODO this could be done by the destructors
    colvars->write_traj_label (colvars->cv_traj_os);
    return COLVARSCRIPT_OK;
  }

  if (argc >= 4) {
    std::string param = argv[3];

    result = "Syntax error";
    return COLVARSCRIPT_ERROR;
  }
  result = "Syntax error";
  return COLVARSCRIPT_ERROR;
}
