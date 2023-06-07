#include <iostream>
#include <string>

#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvarproxy.h"

#include "colvarproxy_stub.h"


extern "C" int main(int argc, char *argv[]) {
  if (argc > 1) {
    colvarproxy_stub *proxy = new colvarproxy_stub();
    int error_code = COLVARS_OK;

    // Initialize simple unit system to test file input
    error_code |= proxy->set_unit_system("real", false);
    error_code |= proxy->colvars->read_config_file(argv[1]);

    std::cout << std::endl
              << "Input files read during this test:"
              << std::endl;
    unsigned char * args[2] = {
      (unsigned char *) "cv",
      (unsigned char *) "listinputfiles" };
    error_code |= run_colvarscript_command(2, args);
    std::cout << get_colvarscript_result() << std::endl;

  } else {
    std::cerr << "ERROR: Missing configuration file." << std::endl;
  }
}
