#include <iostream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"


extern "C" int main(int argc, char *argv[]) {

  colvarproxy *proxy = new colvarproxy();
  proxy->colvars = new colvarmodule(proxy);

  proxy->tcl_run_script("puts \"\n(Tcl) Tcl script running successfully using embedded interpreter.\"");
  proxy->tcl_run_script("puts \"(Tcl) Currently defined variables:\n(Tcl) [info vars]\"");
  proxy->tcl_run_script("puts \"(Tcl) Currently defined procedures:\n(Tcl) [info procs]\"");

  return 0;
}
