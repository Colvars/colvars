#include <iostream>
#include <string>

#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvarproxy.h"

#include "colvarproxy_stub.h"


extern "C" int main(int argc, char *argv[]) {
  if (argc > 1) {
    colvarproxy_stub *proxy = new colvarproxy_stub();
    // Initialize simple unit system to test file input
    proxy->set_unit_system("real", false);
    return proxy->colvars->read_config_file(argv[1]);
  } else {
    std::cerr << "ERROR: Missing configuration file." << std::endl;
  }
}
