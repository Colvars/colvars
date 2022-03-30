#include <iostream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"


extern "C" int main(int argc, char *argv[]) {

  colvarproxy *proxy = new colvarproxy();
  proxy->colvars = new colvarmodule(proxy);

  int res = proxy->tcl_run_script("puts \"\n(Tcl) Tcl script running successfully using embedded interpreter.\"");

  if (res != COLVARS_OK) {
    std::cout << "Error running Tcl script.\n";
    if (res == COLVARS_NOT_IMPLEMENTED)
      std::cout << "Tcl scripts are not enabled in this build of the Colvars library.\n";
    return 1;
  }

  proxy->tcl_run_script("puts \"(Tcl) Currently defined variables:\n(Tcl) [info vars]\"");
  proxy->tcl_run_script("puts \"(Tcl) Currently defined procedures:\n(Tcl) [info procs]\"");

  return 0;
}
