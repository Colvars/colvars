// #define COLVARS_DEBUG true
#include "colvar_gpu_support.h"
#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvarproxy.h"
#include "colvarproxy_stub_gpu.h"
#include "CLI11.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>


int main(int argc, char *argv[]) {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
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

  int err = COLVARS_OK;
  colvarproxy_stub_gpu *proxy = new colvarproxy_stub_gpu();
  proxy->init_cvm();
  colvarmodule* const cvmodule = proxy->cvmodule;
  // Initialize simple unit system to test file input
  err |= proxy->set_unit_system("real", false);
  // Initialize simple unit system to test file input
  err |= proxy->set_unit_system("real", false);

  if (argc > 3) {
    err |= proxy->set_output_prefix(output_prefix);
  }
  err |= proxy->cvmodule->setup_input();
  err |= proxy->cvmodule->setup_output();
  err |= proxy->cvmodule->read_config_file(configuration_file.c_str());
  if (err != COLVARS_OK) {
    cvmodule->log("Error occurred!\n");
  }

  if (trajectory_file.size() > 0) {
    int io_err = 0;
    while (!io_err) {
      io_err = proxy->read_frame_xyz(trajectory_file.c_str(), output_force);
      err = cvmodule->get_error();
      if (err != COLVARS_OK) {
        cvmodule->log("Error occurred!\n");
      }
      if (!io_err) cvmodule->log("Frame " + cvm::to_str(cvmodule->step_absolute()));
    }
    proxy->post_run();
    cvmodule->log("Done");
    err = cvmodule->get_error();
    if (err != COLVARS_OK) {
      cvmodule->log("Error occurred!\n");
    }
  }

  cvmodule->log("Input files read during this test:");
  unsigned char * args[2] = {
    (unsigned char *) "cv",
    (unsigned char *) "listinputfiles" };
  err |= run_colvarscript_command(proxy->script, 2, args);
  cvmodule->log("  " + std::string(get_colvarscript_result()));

  double const max_gradient_error = proxy->cvmodule->get_max_gradient_error();
  if (max_gradient_error > 0.) {
    cvmodule->log("Max gradient error (debugGradients): " + cvm::to_str(max_gradient_error));

    double threshold = 1e-3;
    // Fail test if error is above threshold
    if (max_gradient_error > threshold) {
      cvmodule->log("Error: gradient inaccuracy is above threshold (" + cvm::to_str(threshold) + ")");
      err = 1;
    }
  }

  delete proxy;
  return err;
#else
  std::cout << "This program requires CUDA to test." << std::endl;
  return 1;
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
}
