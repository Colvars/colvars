// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <cstdlib>
#include <cstring>
#include <sstream>

#if defined(NAMD_TCL) || defined(VMDTCL)
#define COLVARS_TCL
#include <tcl.h>
#endif

#include "colvarproxy.h"
#include "colvardeps.h"
#include "colvarscript.h"
#include "colvarscript_commands.h"



colvarscript::colvarscript(colvarproxy *p)
 : proxy_(p),
   colvars(p->colvars),
   proxy_error(0)
{
  cmd_names = NULL;
  init_commands();
#ifdef COLVARS_TCL
  // TODO put this in backend functions so we don't have to delete
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(proxy_->get_tcl_interp());
  Tcl_DeleteCommand(interp, "cv");
  Tcl_CreateObjCommand(interp, "cv", tcl_run_colvarscript_command,
                       (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  cvm::log("Redefining the Tcl \"cv\" command to the new script interface.");
#endif
}


colvarscript::~colvarscript()
{
  if (cmd_names) {
    delete [] cmd_names;
    cmd_names = NULL;
  }
}


int colvarscript::init_commands()
{
  if (cvm::debug()) {
    cvm::log("Called colvarcript::init_commands()\n");
  }

  cmd_help.resize(colvarscript::cv_n_commands);
  cmd_n_args_min.resize(colvarscript::cv_n_commands);
  cmd_n_args_max.resize(colvarscript::cv_n_commands);
  cmd_arghelp.resize(colvarscript::cv_n_commands);
  cmd_fns.resize(colvarscript::cv_n_commands);

  if (cmd_names) {
    delete [] cmd_names;
    cmd_names = NULL;
  }
  cmd_names = new char const * [colvarscript::cv_n_commands];

#undef COLVARSCRIPT_COMMANDS_H // disable include guard
#if defined(CVSCRIPT)
#undef CVSCRIPT // disable default macro
#endif
#define CVSCRIPT_COMM_INIT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGHELP) {   \
    init_command(COMM,#COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGHELP,&(CVSCRIPT_COMM_FNAME(COMM))); \
  }
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_INIT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS)

#include "colvarscript_commands.h"

#undef CVSCRIPT_COMM_INIT
#undef CVSCRIPT

  return COLVARS_OK;
}


int colvarscript::init_command(colvarscript::command const &comm,
                               char const *name, char const *help,
                               int n_args_min, int n_args_max,
                               char const *arghelp,
                               int (*fn)(void *, int, unsigned char * const *))
{
  cmd_str_map[std::string(name)] = comm;
  cmd_names[comm] = name;
  cmd_help[comm] = help;
  cmd_n_args_min[comm] = n_args_min;
  cmd_n_args_max[comm] = n_args_max;
  std::string const arghelp_str(arghelp);
  std::istringstream is(arghelp_str);
  std::string line;
  for (int iarg = 0; iarg < n_args_max; iarg++) {
    if (! std::getline(is, line)) {
      return cvm::error("Error: could not initialize help string for scripting "
                        "command \""+std::string(name)+"\".\n", BUG_ERROR);
    }
    cmd_arghelp[comm].push_back(line);
  }
  cmd_fns[comm] = fn;
  if (cvm::debug()) {
    cvm::log("Defined command \""+std::string(name)+"\", with help string:\n");
    cvm::log(get_command_help(name));
  }
  return COLVARS_OK;
}


std::string colvarscript::get_command_help(char const *cmd)
{
  if (cmd_str_map.count(cmd) > 0) {
    colvarscript::command const c = cmd_str_map[std::string(cmd)];
    std::string new_result(cmd_help[c]+"\n");
    if (cmd_n_args_max[c] == 0) return new_result;
    new_result += "\nParameters\n";
    new_result += "----------\n\n";
    size_t i;
    for (i = 0; i < cmd_n_args_min[c]; i++) {
      new_result += cmd_arghelp[c][i]+"\n\n";
    }
    for (i = cmd_n_args_min[c]; i < cmd_n_args_max[c]; i++) {
      new_result += cmd_arghelp[c][i]+" (optional)\n\n";
    }
    return new_result;
  }

  cvm::error("Error: command "+std::string(cmd)+
             " is not implemented.\n", INPUT_ERROR);
  return std::string("");
}


int colvarscript::run(int objc, unsigned char *const objv[])
{
  result.clear();

  if (cvm::debug()) {
    cvm::log("Called script run with " + cvm::to_str(objc) + " args:");
    for (int i = 0; i < objc; i++) {
      cvm::log(obj_to_str(objv[i]));
    }
  }

  if (objc < 2) {
    set_result_str("No commands given: use \"cv help\" "
                   "for a list of commands.");
    return COLVARSCRIPT_ERROR;
  }

  // Name of the (sub)command
  std::string const cmd(obj_to_str(objv[1]));

  // Build a safe-to-print command line to print in case of error
  std::string cmdline(std::string(obj_to_str(objv[0]))+std::string(" ")+cmd);

  // Pointer to the function implementing it
  int (*cmd_fn)(void *, int, unsigned char * const *) = NULL;

  // Pointer to object handling the command (the functions are C)
  void *obj_for_cmd = NULL;

  if (cmd == "colvar") {

    if (objc < 4) {
      add_error_msg("Missing parameters\n" + help_string());
      return COLVARSCRIPT_ERROR;
    }
    std::string const name(obj_to_str(objv[2]));
    obj_for_cmd = reinterpret_cast<void *>(cvm::colvar_by_name(name));
    if (obj_for_cmd == NULL) {
      add_error_msg("Colvar not found: " + name);
      return COLVARSCRIPT_ERROR;
    }
    std::string const &subcmd(obj_to_str(objv[3]));
    cmd_fn = get_cmd_fn(std::string("colvar_")+subcmd);
    cmdline += std::string(" <name> ")+subcmd;
    if (objc > 4) {
      cmdline += " ...";
    }

  } else if (cmd == "bias") {

    if (objc < 4) {
      add_error_msg("Missing parameters\n" + help_string());
      return COLVARSCRIPT_ERROR;
    }
    std::string const name(obj_to_str(objv[2]));
    obj_for_cmd = reinterpret_cast<void *>(cvm::bias_by_name(name));
    if (obj_for_cmd == NULL) {
      add_error_msg("Bias not found: " + name);
      return COLVARSCRIPT_ERROR;
    }
    std::string const &subcmd(obj_to_str(objv[3]));
    cmd_fn = get_cmd_fn(std::string("bias_")+subcmd);
    cmdline += std::string(" <name> ")+subcmd;
    if (objc > 4) {
      cmdline += " ...";
    }

  } else {

    cmd_fn = get_cmd_fn(std::string(std::string("cv_"+cmd)));
    obj_for_cmd = reinterpret_cast<void *>(this);

    if (objc > 2) {
      cmdline += " ...";
    }
  }

  int error_code = COLVARS_OK;

  // If command was found in map, execute it
  if (cmd_fn) {
    error_code = (*cmd_fn)(obj_for_cmd, objc, objv);
  } else {
    add_error_msg("Syntax error: "+cmdline+"\n"
                  "  Run \"cv listcommands\" or \"cv help <command>\" "
                  "to get the correct syntax.\n");
    error_code = COLVARSCRIPT_ERROR;
  }

  return error_code;
}


int colvarscript::proc_features(colvardeps *obj,
                                int objc, unsigned char *const objv[]) {

  // size was already checked before calling
  std::string const subcmd(obj_to_str(objv[3]));

  if (cvm::debug()) {
    cvm::log("Called proc_features() with " + cvm::to_str(objc) + " args:");
    for (int i = 0; i < objc; i++) {
      cvm::log(obj_to_str(objv[i]));
    }
  }

  if ((subcmd == "get") || (subcmd == "set")) {
    std::vector<colvardeps::feature *> const &features = obj->features();
    std::string const req_feature(obj_to_str(objv[4]));
    colvardeps::feature *f = NULL;
    int fid = 0;
    for (fid = 0; fid < int(features.size()); fid++) {
      if (features[fid]->description ==
          colvarparse::to_lower_cppstr(req_feature)) {
        f = features[fid];
        break;
      }
    }

    if (f == NULL) {

      add_error_msg("Error: feature \""+req_feature+"\" does not exist.\n");
      return COLVARSCRIPT_ERROR;

    } else {

      if (! obj->is_available(fid)) {
        add_error_msg("Error: feature \""+req_feature+"\" is unavailable.\n");
        return COLVARSCRIPT_ERROR;
      }

      if (subcmd == "get") {
        set_result_str(cvm::to_str(obj->is_enabled(fid) ? 1 : 0));
        return COLVARS_OK;
      }

      if (subcmd == "set") {
        if (objc == 6) {
          std::string const yesno =
            colvarparse::to_lower_cppstr(std::string(obj_to_str(objv[5])));
          if ((yesno == std::string("yes")) ||
              (yesno == std::string("on")) ||
              (yesno == std::string("1"))) {
            obj->enable(fid);
            return COLVARS_OK;
          } else if ((yesno == std::string("no")) ||
              (yesno == std::string("off")) ||
              (yesno == std::string("0"))) {
            obj->disable(fid);
            return COLVARS_OK;
          }
        }
        add_error_msg("Missing value when setting feature \""+req_feature+
                      "\".\n");
        return COLVARSCRIPT_ERROR;
      }
    }
  }

  // This shouldn't be reached any more
  return COLVARSCRIPT_ERROR;
}


std::string colvarscript::help_string() const
{
  std::string buf;
  buf = "Usage: cv <subcommand> [args...]\n\
\n\
Managing the Colvars module:\n\
  configfile <file name>      -- read configuration from a file\n\
  config <string>             -- read configuration from the given string\n\
  getconfig                   -- get the module's configuration string\n\
  resetindexgroups            -- clear the index groups loaded so far\n\
  reset                       -- delete all internal configuration\n\
  delete                      -- delete this Colvars module instance\n\
  version                     -- return version of Colvars Module\n\
  \n\
Input and output:\n\
  list                        -- return a list of all variables\n\
  list biases                 -- return a list of all biases\n\
  load <file name>            -- load a state file (requires configuration)\n\
  save <file name>            -- save a state file (requires configuration)\n\
  update                      -- recalculate colvars and biases\n\
  addenergy <E>               -- add <E> to the total bias energy\n\
  printframe                  -- return a summary of the current frame\n\
  printframelabels            -- return labels to annotate printframe's output\n";

  long int tmp;
  if (proxy_->get_frame(tmp) != COLVARS_NOT_IMPLEMENTED) {
      buf += "\
  frame                       -- return current frame number\n\
  frame <new_frame>           -- set frame number\n";
  }

  buf += "\n\
Accessing collective variables:\n\
  colvar <name> value         -- return the current value of colvar <name>\n\
  colvar <name> update        -- recalculate colvar <name>\n\
  colvar <name> type          -- return the type of colvar <name>\n\
  colvar <name> delete        -- delete colvar <name>\n\
  colvar <name> addforce <F>  -- apply given force on colvar <name>\n\
  colvar <name> getappliedforce -- return applied force of colvar <name>\n\
  colvar <name> gettotalforce -- return total force of colvar <name>\n\
  colvar <name> getconfig     -- return config string of colvar <name>\n\
  colvar <name> cvcflags <fl> -- enable or disable cvcs according to 0/1 flags\n\
  colvar <name> modifycvcs <str> -- pass new config strings to each CVC\n\
  colvar <name> get <f>       -- get the value of the colvar feature <f>\n\
  colvar <name> set <f> <val> -- set the value of the colvar feature <f>\n\
\n\
Accessing biases:\n\
  bias <name> energy          -- return the current energy of bias <name>\n\
  bias <name> update          -- recalculate bias <name>\n\
  bias <name> delete          -- delete bias <name>\n\
  bias <name> getconfig       -- return config string of bias <name>\n\
  bias <name> get <f>         -- get the value of the bias feature <f>\n\
  bias <name> set <f> <val>   -- set the value of the bias feature <f>\n\
";

  return buf;
}


int colvarscript::unsupported_op()
{
  return cvm::error("Error: unsupported script operation.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarscript::set_result_str(std::string const &s)
{
  result = s;
  return COLVARS_OK;
}


void colvarscript::add_error_msg(std::string const &s)
{
  result += s;
  // Ensure terminating newlines
  if (s[s.size()-1] != '\n') {
    result += "\n";
  }
}


int colvarscript::clear_str_result()
{
  result.clear();
  return COLVARS_OK;
}


extern "C"
int run_colvarscript_command(int objc, unsigned char *const objv[])
{
  colvarmodule *cv = cvm::main();
  colvarscript *script = cv ? cv->proxy->script : NULL;
  if (!script) {
    cvm::error("Called run_colvarscript_command without a script object.\n",
               BUG_ERROR);
    return -1;
  }
  int retval = script->run(objc, objv);
  return retval;
}


extern "C"
const char * get_colvarscript_result()
{
  colvarscript *script = colvarscript_obj();
  if (!script) {
    cvm::error("Called get_colvarscript_result without a script object.\n");
    return NULL;
  }
  return script->str_result().c_str();
}


#if defined(COLVARS_TCL)

/// Run the script API via Tcl command-line interface
/// \param clientData Not used
/// \param my_interp Pointer to Tcl_Interp object (read from Colvars if NULL)
/// \param objc Number of Tcl command parameters
/// \param objv Array of command parameters
/// \return Result of the script command
extern "C"
int tcl_run_colvarscript_command(ClientData clientData,
                                 Tcl_Interp *my_interp,
                                 int objc, Tcl_Obj *const objv[])
{
  colvarmodule *colvars = cvm::main();

  if (!colvars) {
    Tcl_SetResult(my_interp, const_cast<char *>("Colvars module not active"),
                  TCL_VOLATILE);
    return TCL_ERROR;
  }

  colvarproxy *proxy = colvars->proxy;
  Tcl_Interp *interp = my_interp ? my_interp :
    reinterpret_cast<Tcl_Interp *>(proxy->get_tcl_interp());
  colvarscript *script = colvarscript_obj();
  if (!script) {
    char const *errstr = "Called tcl_run_colvarscript_command "
      "without a Colvars script interface set up.\n";
    Tcl_SetResult(interp, const_cast<char *>(errstr), TCL_VOLATILE);
    return TCL_ERROR;
  }

  cvm::clear_error();

  int retval = script->run(objc,
                           reinterpret_cast<unsigned char * const *>(objv));

  std::string result = proxy->get_error_msgs() + script->result;

  Tcl_SetResult(interp, const_cast<char *>(result.c_str()),
                TCL_VOLATILE);

  return (retval == COLVARS_OK) ? TCL_OK : TCL_ERROR;
}

#endif // #if defined(COLVARS_TCL)
