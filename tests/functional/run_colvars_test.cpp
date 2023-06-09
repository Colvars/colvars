#include <iostream>
#include <fstream>
#include <string>

#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"

extern "C" int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 4) {
    std::cerr << "Wrong number of arguments.\n"
              << "Usage: run_colvars_test <configuration_file> [XYZ_trajectory_file] [output_prefix]"
              << std::endl;
    return 1;
  }
  int err = 0;

  colvarproxy_stub *proxy = new colvarproxy_stub();
  // Initialize simple unit system to test file input
  err |= proxy->set_unit_system("real", false);

  if (argc > 3) {
    err |= proxy->set_output_prefix(argv[3]);
  }
  err |= proxy->colvars->setup_input();
  err |= proxy->colvars->setup_output();
  err |= proxy->colvars->read_config_file(argv[1]);

  if (argc > 2) {
    // Read number of atoms from XYZ header
    std::ifstream ifs(argv[2]);
    int natoms;
    ifs >> natoms;
    ifs.close();
    cvm::log("Reading trajectory for " + cvm::to_str(natoms)
              + " atoms from XYZ file " + argv[2]);
    for (int ai = 0; ai < natoms; ai++) {
      err |= proxy->init_atom(ai+1);
    }
    int io_err = 0;
    while (!io_err) {
      io_err = proxy->read_frame_xyz(argv[2]);
      if (!io_err) cvm::log("Frame " + cvm::to_str(cvm::step_absolute()));
    }
    proxy->post_run();
    cvm::log("Done");
  }

  cvm::log("Input files read during this test:");
  unsigned char * args[2] = {
    (unsigned char *) "cv",
    (unsigned char *) "listinputfiles" };
  err |= run_colvarscript_command(2, args);
  cvm::log("  " + std::string(get_colvarscript_result()));

  return err;
}
