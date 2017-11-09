// -*- c++ -*-

#include <cstdlib>
#include <stdlib.h>
#include <string.h>

#if defined(NAMD_TCL) || defined(VMDTCL)
#define COLVARS_TCL
#include <tcl.h>
#endif

#include "colvarproxy.h"
#include "colvardeps.h"
#include "colvarscript.h"
#include "colvarscript_commands.h"



colvarscript::colvarscript(colvarproxy *p)
 : proxy(p),
   colvars(p->colvars),
   proxy_error(0)
{
  comm_names = NULL;
  init_commands();
}


colvarscript::~colvarscript()
{
  if (comm_names) {
    delete [] comm_names;
    comm_names = NULL;
  }
}


int colvarscript::init_commands()
{
  comm_help.resize(colvarscript::cv_n_commands);
  comm_n_args_min.resize(colvarscript::cv_n_commands);
  comm_n_args_max.resize(colvarscript::cv_n_commands);
  comm_arghelp.resize(colvarscript::cv_n_commands);
  comm_fns.resize(colvarscript::cv_n_commands);

  if (comm_names) {
    delete [] comm_names;
    comm_names = NULL;
  }
  comm_names = new char const * [colvarscript::cv_n_commands];

#undef COLVARSCRIPT_COMMANDS_H // disable include guard
#if defined(CVSCRIPT)
#undef CVSCRIPT // disable default macro
#endif
#define CVSCRIPT_COMM_INIT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS) {      \
    comm_str_map[#COMM] = COMM;                                         \
    comm_names[COMM] = #COMM;                                           \
    comm_help[COMM] = HELP;                                             \
    comm_n_args_min[COMM] = N_ARGS_MIN;                                 \
    comm_n_args_max[COMM] = N_ARGS_MAX;                                 \
    comm_arghelp[COMM] = ARGS;                                          \
    comm_fns[COMM] = &(CVSCRIPT_COMM_FNAME(COMM));                      \
  }
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_INIT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS)

#include "colvarscript_commands.h"

#undef CVSCRIPT_COMM_INIT
#undef CVSCRIPT

  return COLVARS_OK;
}


std::string colvarscript::get_command_help(char const *cmd)
{
  if (comm_str_map.count(cmd) > 0) {
    colvarscript::command const c = comm_str_map[std::string(cmd)];
    std::string result(comm_help[c]+"\n");
    for (size_t i = 0; i < comm_n_args_max[c]; i++) {
      result += comm_arghelp[c][i]+"\n";
    }
    return result;
  }

  cvm::error("Error: command "+std::string(cmd)+
             " is not implemented.\n", INPUT_ERROR);
  return std::string("");
}


// int (*colvarscript::get_command(colvarscript::command c))(
//    void *, int, unsigned char * const *)
// {
//   return comm_fns[c];
// }


int colvarscript::run(int objc, unsigned char *const objv[])
{
  colvarscript *script = this;
  colvarproxy *proxy = cvm::main()->proxy;
#if defined(COLVARS_TCL)
  Tcl_Interp *interp =
    reinterpret_cast<Tcl_Interp *>(proxy->get_tcl_interp());
#endif

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

  std::string const cmd(obj_to_str(objv[1]));

  int error_code = COLVARS_OK;

  if (cmd == "newstyle") {
#if defined(COLVARS_TCL)
    Tcl_DeleteCommand(interp, "cv");
    Tcl_CreateObjCommand(interp, "cv", tcl_run_colvarscript_command,
                         (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
    cvm::log("Switching to new scripting interface.");
#endif
    return COLVARS_OK;
  }

  // If new-style command is found in map, execute it
  std::string const cmd_key("cv_"+cmd);
  if (comm_str_map.count(cmd_key) > 0) {
    error_code |= (*(comm_fns[comm_str_map[cmd_key]]))(
                      reinterpret_cast<void *>(this), objc, objv);
    return error_code;
  }

  if (cmd == "colvar") {
    if (objc < 3) {
      result = "Missing parameters\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    std::string const name(obj_to_str(objv[2]));
    colvar *cv = cvm::colvar_by_name(name);
    if (cv == NULL) {
      result = "Colvar not found: " + name;
      return COLVARSCRIPT_ERROR;
    }
    return proc_colvar(cv, objc-1, &(objv[1]));
  }

  if (cmd == "bias") {
    if (objc < 3) {
      result = "Missing parameters\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    std::string const name(obj_to_str(objv[2]));
    colvarbias *b = cvm::bias_by_name(name);
    if (b == NULL) {
      result = "Bias not found: " + name;
      return COLVARSCRIPT_ERROR;
    }
    return proc_bias(b, objc-1, &(objv[1]));
  }

  if (cmd == "version") {
    result = COLVARS_VERSION;
    return COLVARS_OK;
  }

  if (cmd == "reset") {
    /// Delete every child object
    colvars->reset();
    return COLVARS_OK;
  }

  if (cmd == "delete") {
    // Note: the delete bit may be ignored by some backends
    // it is mostly useful in VMD
    return proxy->request_deletion();
  }

  if (cmd == "update") {
    error_code |= proxy->update_input();
    error_code |= colvars->calc();
    error_code |= proxy->update_output();
    if (error_code) {
      result += "Error updating the Colvars module.\n";
    }
    return error_code;
  }

  if (cmd == "list") {
    if (objc == 2) {
      for (std::vector<colvar *>::iterator cvi = colvars->colvars.begin();
           cvi != colvars->colvars.end();
           ++cvi) {
        result += (cvi == colvars->colvars.begin() ? "" : " ") + (*cvi)->name;
      }
      return COLVARS_OK;
    } else if (objc == 3 && !strcmp(obj_to_str(objv[2]), "biases")) {
      for (std::vector<colvarbias *>::iterator bi = colvars->biases.begin();
           bi != colvars->biases.end();
           ++bi) {
        result += (bi == colvars->biases.begin() ? "" : " ") + (*bi)->name;
      }
      return COLVARS_OK;
    } else {
      result = "Wrong arguments to command \"list\"\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Parse config from file
  if (cmd == "configfile") {
    if (objc < 3) {
      result = "Missing arguments\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    if (colvars->read_config_file(obj_to_str(objv[2])) == COLVARS_OK) {
      return COLVARS_OK;
    } else {
      result = "Error parsing configuration file";
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Parse config from string
  if (cmd == "config") {
    return exec_command(cv_config, NULL, objc, objv);
  }

  /// Load an input state file
  if (cmd == "load") {
    if (objc < 3) {
      result = "Missing arguments\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    proxy->input_prefix() = obj_to_str(objv[2]);
    if (colvars->setup_input() == COLVARS_OK) {
      return COLVARS_OK;
    } else {
      result = "Error loading state file";
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Save to an output state file
  if (cmd == "save") {
    if (objc < 3) {
      result = "Missing arguments";
      return COLVARSCRIPT_ERROR;
    }
    proxy->output_prefix() = obj_to_str(objv[2]);
    int error = 0;
    error |= colvars->setup_output();
    error |= colvars->write_restart_file(colvars->output_prefix()+
                                         ".colvars.state");
    error |= colvars->write_output_files();
    return error ? COLVARSCRIPT_ERROR : COLVARS_OK;
  }

  /// Print the values that would go on colvars.traj
  if (cmd == "printframelabels") {
    std::ostringstream os;
    colvars->write_traj_label(os);
    result = os.str();
    return COLVARS_OK;
  }
  if (cmd == "printframe") {
    std::ostringstream os;
    colvars->write_traj(os);
    result = os.str();
    return COLVARS_OK;
  }

  if (cmd == "frame") {
    if (objc == 2) {
      long int f;
      int error = proxy->get_frame(f);
      if (error == COLVARS_OK) {
        result = cvm::to_str(f);
        return COLVARS_OK;
      } else {
        result = "Frame number is not available";
        return COLVARSCRIPT_ERROR;
      }
    } else if (objc == 3) {
      // Failure of this function does not trigger an error, but
      // returns nonzero, to let scripts detect available frames
      int error = proxy->set_frame(strtol(obj_to_str(objv[2]), NULL, 10));
      result = cvm::to_str(error == COLVARS_OK ? 0 : -1);
      return COLVARS_OK;
    } else {
      result = "Wrong arguments to command \"frame\"\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
  }

  if (cmd == "addenergy") {
    if (objc == 3) {
      colvars->total_bias_energy += strtod(obj_to_str(objv[2]), NULL);
      return COLVARS_OK;
    } else {
      result = "Wrong arguments to command \"addenergy\"\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
  }

  if (cmd == "help") {
    return exec_command(cv_help, NULL, objc, objv);
  }

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_colvar(colvar *cv, int objc, unsigned char *const objv[]) {

  std::string const subcmd(obj_to_str(objv[2]));

  if (subcmd == "value") {
    result = (cv->value()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "run_ave") {
    result = (cv->run_ave()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "width") {
    result = cvm::to_str(cv->width, 0, cvm::cv_prec);
    return COLVARS_OK;
  }

  if (subcmd == "type") {
    result = cv->value().type_desc(cv->value().value_type);
    return COLVARS_OK;
  }

  if (subcmd == "update") {
    cv->calc();
    cv->update_forces_energy();
    result = (cv->value()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "delete") {
    size_t i;
    for (i = 0; i < cv->biases.size(); i++) {
      delete cv->biases[i];
    }
    cv->biases.clear();
    // colvar destructor is tasked with the cleanup
    delete cv;
    // TODO this could be done by the destructors
    if (colvars->cv_traj_os != NULL) {
      colvars->write_traj_label(*(colvars->cv_traj_os));
    }
    return COLVARS_OK;
  }

  if (subcmd == "getconfig") {
    result = cv->get_config();
    return COLVARS_OK;
  }

  if (subcmd == "getappliedforce") {
    result = (cv->applied_force()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "getsystemforce") {
    // TODO warning here
    result = (cv->total_force()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "gettotalforce") {
    result = (cv->total_force()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "addforce") {
    if (objc < 4) {
      result = "addforce: missing parameter: force value\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    std::string const f_str(obj_to_str(objv[3]));
    std::istringstream is(f_str);
    is.width(cvm::cv_width);
    is.precision(cvm::cv_prec);
    colvarvalue force(cv->value());
    force.is_derivative();
    if (force.from_simple_string(is.str()) != COLVARS_OK) {
      result = "addforce : error parsing force value";
      return COLVARSCRIPT_ERROR;
    }
    cv->add_bias_force(force);
    result = force.to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "cvcflags") {
    if (objc < 4) {
      result = "cvcflags: missing parameter: vector of flags";
      return COLVARSCRIPT_ERROR;
    }
    std::string const flags_str(obj_to_str(objv[3]));
    std::istringstream is(flags_str);
    std::vector<bool> flags;

    int flag;
    while (is >> flag) {
      flags.push_back(flag != 0);
    }

    int res = cv->set_cvc_flags(flags);
    if (res != COLVARS_OK) {
      result = "Error setting CVC flags";
      return COLVARSCRIPT_ERROR;
    }
    result = "0";
    return COLVARS_OK;
  }

  if ((subcmd == "get") || (subcmd == "set") || (subcmd == "state")) {
    return proc_features(cv, objc, objv);
  }

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_bias(colvarbias *b, int objc, unsigned char *const objv[]) {

  std::string const subcmd(obj_to_str(objv[2]));

  if (subcmd == "energy") {
    result = cvm::to_str(b->get_energy());
    return COLVARS_OK;
  }

  if (subcmd == "update") {
    b->update();
    result = cvm::to_str(b->get_energy());
    return COLVARS_OK;
  }

  if (subcmd == "getconfig") {
    result = b->get_config();
    return COLVARS_OK;
  }

  // Subcommands for MW ABF
  if (subcmd == "bin") {
    int r = b->current_bin();
    result = cvm::to_str(r);
    return COLVARS_OK;
  }

  if (subcmd == "binnum") {
    int r = b->bin_num();
    if (r < 0) {
      result = "Error: calling bin_num() for bias " + b->name;
      return COLVARSCRIPT_ERROR;
    }
    result = cvm::to_str(r);
    return COLVARS_OK;
  }

  if (subcmd == "share") {
    int r = b->replica_share();
    if (r < 0) {
      result = "Error: calling replica_share() for bias " + b->name;
      return COLVARSCRIPT_ERROR;
    }
    result = cvm::to_str(r);
    return COLVARS_OK;
  }
  // End commands for MW ABF

  if (subcmd == "delete") {
    // the bias destructor takes care of the cleanup at cvm level
    delete b;
    // TODO this could be done by the destructors
    if (colvars->cv_traj_os != NULL) {
      colvars->write_traj_label(*(colvars->cv_traj_os));
    }
    return COLVARS_OK;
  }

  if ((subcmd == "get") || (subcmd == "set") || (subcmd == "state")) {
    return proc_features(b, objc, objv);
  }

  if (objc >= 4) {
    std::string const param(obj_to_str(objv[3]));
    if (subcmd == "count") {
      int index;
      if (!(std::istringstream(param) >> index)) {
        result = "bin_count: error parsing bin index";
        return COLVARSCRIPT_ERROR;
      }
      result = cvm::to_str(b->bin_count(index));
      return COLVARS_OK;
    }

    result = "Syntax error\n" + help_string();
    return COLVARSCRIPT_ERROR;
  }

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_features(colvardeps *obj,
                                int objc, unsigned char *const objv[]) {
  // size was already checked before calling
  std::string const subcmd(obj_to_str(objv[2]));

  if (objc == 3) {
    if (subcmd == "state") {
      // TODO make this returned as result?
      obj->print_state();
      return COLVARS_OK;
    }

    // get and set commands require more arguments
    result = "Syntax error\n" + help_string();
    return COLVARSCRIPT_ERROR;
  }

  if ((subcmd == "get") || (subcmd == "set")) {
    std::vector<colvardeps::feature *> const &features = obj->features();
    std::string const req_feature(obj_to_str(objv[3]));
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

      result = "Error: feature \""+req_feature+"\" does not exist.\n";
      return COLVARSCRIPT_ERROR;

    } else {

      if (! obj->is_available(fid)) {
        result = "Error: feature \""+req_feature+"\" is unavailable.\n";
        return COLVARSCRIPT_ERROR;
      }

      if (subcmd == "get") {
        result = cvm::to_str(obj->is_enabled(fid) ? 1 : 0);
        return COLVARS_OK;
      }

      if (subcmd == "set") {
        if (objc == 5) {
          std::string const yesno =
            colvarparse::to_lower_cppstr(std::string(obj_to_str(objv[4])));
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
        result = "Syntax error\n" + help_string();
        return COLVARSCRIPT_ERROR;
      }
    }
  }

  result = "Syntax error\n" + help_string();
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
  reset                       -- delete all internal configuration\n\
  delete                      -- delete this Colvars module instance\n\
  version                     -- return version of Colvars code\n\
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
  if (proxy->get_frame(tmp) != COLVARS_NOT_IMPLEMENTED) {
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
                    COLVARS_ERROR); // TODO make an entry for script
}


int colvarscript::set_result_str(std::string const &s)
{
  result = s;
  return COLVARS_OK;
}


int colvarscript::clear_str_result()
{
  result.clear();
  return COLVARS_OK;
}


int colvarscript::set_result_real(cvm::real x, unsigned char *obj)
{
  colvarproxy *proxy = cvm::main()->proxy;
  switch (interp_type_list.back()) {
  case cv_text:
    {
      std::string const x_str = cvm::to_str(x, cvm::cv_width, cvm::cv_prec);
      if (obj) {
        char *obj_str = reinterpret_cast<char *>(obj);
        strcpy(obj_str, x_str.c_str());
      } else {
        set_result_str(x_str);
      }
    }
    break;
  case cv_nointerp:
    if (obj) {
      double *y = reinterpret_cast<double *>(obj);
      *y = x;
    } else {
      return unsupported_op(); // TODO define memory management
    }
    break;
  case cv_tcl:
    // TODO handle obj
    *(modify_obj_result()) =
      reinterpret_cast<unsigned char *>(proxy->tcl_set_real(x));
    *(modify_obj_result_size()) = -1;
    break;
  default:
    return unsupported_op();
    break;
  }
  return COLVARS_OK;
}


int colvarscript::set_result_real_vec(std::vector<cvm::real> const &x,
                                      unsigned char *obj)
{
  colvarproxy *proxy = cvm::main()->proxy;
  switch (interp_type_list.back()) {
  case cv_text:
    {
      std::string const x_str = cvm::to_str(x, cvm::cv_width, cvm::cv_prec);
      if (obj) {
        char *obj_str = reinterpret_cast<char *>(obj);
        strcpy(obj_str, x_str.c_str());
      } else {
        set_result_str(x_str);
      }
    }
    break;
  case cv_nointerp:
    if (obj) {
      double *out = reinterpret_cast<double *>(obj);
      for (size_t i = 0; i < x.size(); i++) {
        out[i] = x[i];
      }
    } else {
      return unsupported_op();
    }
    break;
  case cv_tcl:
    if (obj) {
      unsigned char * tcl_obj =
        reinterpret_cast<unsigned char *>(proxy->tcl_list_from_real_vec(x));
      *(modify_obj_result()) = tcl_obj;
      *(modify_obj_result_size()) = -1;
    } else {
      // TODO
      return unsupported_op();
    }
    break;
  default:
    return unsupported_op();
    break;
  }
  return COLVARS_OK;
}


int colvarscript::set_result_int_vec(std::vector<int> const &x,
                                     unsigned char *obj)
{
  colvarproxy *proxy = cvm::main()->proxy;
  switch (interp_type_list.back()) {
  case cv_text:
    {
      std::string const x_str = cvm::to_str(x, cvm::cv_width, cvm::cv_prec);
      if (obj) {
        char *obj_str = reinterpret_cast<char *>(obj);
        strcpy(obj_str, x_str.c_str());
      } else {
        set_result_str(x_str);
      }
    }
    break;
  case cv_nointerp:
    if (obj) {
      int *out = reinterpret_cast<int *>(obj);
      for (size_t i = 0; i < x.size(); i++) {
        out[i] = x[i];
      }
    } else {
      return unsupported_op();
    }
    break;
  case cv_tcl:
    if (obj) {
      unsigned char * tcl_obj =
        reinterpret_cast<unsigned char *>(proxy->tcl_list_from_int_vec(x));
      *(modify_obj_result()) = tcl_obj;
      *(modify_obj_result_size()) = -1;
    } else {
      // TODO
      return unsupported_op();
    }
    break;
  default:
    return unsupported_op();
    break;
  }
  return COLVARS_OK;
}


#if defined(COLVARS_TCL)

extern "C"
int tcl_run_colvarscript_command(ClientData clientData,
                                 Tcl_Interp *interp_in,
                                 int objc, Tcl_Obj *const objv[])
{
  colvarscript *script = colvarscript_obj();
  colvarproxy *proxy = cvm::main()->proxy;
  Tcl_Interp *interp = interp_in ? interp_in :
    reinterpret_cast<Tcl_Interp *>(proxy->get_tcl_interp());
  if (!script) {
    char const *errstr = "Called tcl_run_colvarscript_command "
      "without a Colvars script interface set up.\n";
    Tcl_SetResult(interp, const_cast<char *>(errstr), TCL_STATIC);
    return TCL_ERROR;
  }
  script->enter_interp_call(colvarscript::cv_tcl);
  int retval = script->run(objc,
                           reinterpret_cast<unsigned char * const *>(objv));
  Tcl_Obj *obj = reinterpret_cast<Tcl_Obj *>(script->obj_result());
  std::string const &str_result = script->str_result();
  if (obj == NULL) {
    obj = Tcl_NewStringObj(str_result.c_str(), str_result.length());
  }
  Tcl_SetObjResult(interp, obj);
  script->exit_interp_call();
  return (retval == COLVARS_OK) ? TCL_OK : TCL_ERROR;
}

#endif // #if defined(COLVARS_TCL)
