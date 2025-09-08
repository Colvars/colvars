#include <iostream>
#include <fstream>
#include <string>

#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"

#include "CLI11.hpp"

int main(int argc, char *argv[]) {
  CLI::App app{"Colvars stub interface for testing"};
  argv = app.ensure_utf8(argv);
  std::string configuration_file;
  std::string output_prefix;
  std::string trajectory_file;
  bool output_force = false;
  app.add_option("-c,--configuration_file", configuration_file, "Input Colvars configuration file")->required(true);
  app.add_option("-t,--trajectory_file", trajectory_file, "Input trajectory file")->required(true);
  const auto output_prefix_option =
    app.add_option("-o,--output_prefix", output_prefix, "Output file prefix")->required(false);
  app.add_flag("--force,!--no-force", output_force, "Write the force files")->needs(output_prefix_option);
  CLI11_PARSE(app, argc, argv);

  int err = 0;

  colvarproxy_stub *proxy = new colvarproxy_stub();
  // Initialize simple unit system to test file input
  err |= proxy->set_unit_system("real", false);

  if (argc > 3) {
    err |= proxy->set_output_prefix(output_prefix);
  }
  err |= proxy->colvars->setup_input();
  err |= proxy->colvars->setup_output();
  err |= proxy->colvars->read_config_file(configuration_file.c_str());

  if (argc > 2) {
    // Read number of atoms from XYZ header
    std::ifstream ifs(trajectory_file);
    int natoms;
    ifs >> natoms;
    ifs.close();
    cvm::log("Reading trajectory for " + cvm::to_str(natoms)
              + " atoms from XYZ file " + trajectory_file);
    for (int ai = 0; ai < natoms; ai++) {
      proxy->init_atom(ai+1);
    }
    int io_err = 0;
    while (!io_err) {
      io_err = proxy->read_frame_xyz(trajectory_file.c_str(), output_force);
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

  double const max_gradient_error = proxy->colvars->get_max_gradient_error();
  if (max_gradient_error > 0.) {
    cvm::log("Max gradient error (debugGradients): " + cvm::to_str(max_gradient_error));

    double threshold = 1e-3;
    // Fail test if error is above threshold
    if (max_gradient_error > threshold) {
      cvm::log("Error: gradient inaccuracy is above threshold (" + cvm::to_str(threshold) + ")");
      err = 1;
    }
  }

  delete proxy;

  return err;
}
