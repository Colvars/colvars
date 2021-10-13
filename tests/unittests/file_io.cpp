#include <iostream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"


extern "C" int main(int argc, char *argv[]) {

  colvarproxy *proxy = new colvarproxy();
  proxy->colvars = new colvarmodule(proxy);
  proxy->script = new colvarscript(proxy);

  proxy->backup_file("nonexistent.txt");

  proxy->remove_file("nonexistent.txt");

  // Produce an (almost) empty state file, twice (uses proxy->backup_file())
  proxy->colvars->write_restart_file("test.colvars.state");
  proxy->colvars->write_restart_file("test.colvars.state");

  proxy->backup_file("test.colvars.state.old");

  proxy->remove_file("test.colvars.state");
  proxy->remove_file("test.colvars.state.old");
  proxy->remove_file("test.colvars.state.old.old");


  return 0;
}
