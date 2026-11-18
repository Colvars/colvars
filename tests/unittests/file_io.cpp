#include <iostream>
#include <fstream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"


int main(int argc, char *argv[]) {

  colvarproxy *proxy = new colvarproxy();
  proxy->cvmodule = new colvarmodule(proxy);

  proxy->backup_file("nonexistent.txt");
  // A missing temporary file must not discard a previous valid snapshot.
  std::ofstream("nonexistent.txt.old") << "saved";
  if (proxy->remove_file("nonexistent.txt") != COLVARS_OK ||
      proxy->rename_file("nonexistent.txt", "nonexistent.txt.old") == COLVARS_OK) {
    return 1;
  }
  {
    std::string contents;
    if (!(std::ifstream("nonexistent.txt.old") >> contents) || contents != "saved") return 1;
  }
  proxy->remove_file("nonexistent.txt.old");

  // Produce an (almost) empty state file, twice (uses proxy->backup_file())
  proxy->cvmodule->write_restart_file("test.colvars.state");
  proxy->cvmodule->write_restart_file("test.colvars.state");

  proxy->backup_file("test.colvars.state.old");

  proxy->remove_file("test.colvars.state");
  proxy->remove_file("test.colvars.state.old");
  proxy->remove_file("test.colvars.state.old.old");


  return 0;
}
