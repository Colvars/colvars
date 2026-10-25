#include <iostream>
#include <fstream>
#include <string>

#include "colvarmodule.h"
#include "colvarscript.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"

#include "CLI11.hpp"

#if defined(COLVARS_TRAP_FPE)
#include <cfenv>
// Enable common floating-point exceptions on Linux/glibc
void enable_fp_traps() {
    std::feclearexcept(FE_ALL_EXCEPT);
#ifdef __GLIBC__
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    std::cout << "Floating-point traps enabled (glibc)\n";
#else
    std::cout << "feenableexcept() not available on this platform\n";
#endif
}
#endif


int main(int argc, char *argv[]) {

  #if defined(COLVARS_TRAP_FPE)
  // Trap floating-point exceptions
  enable_fp_traps();
  #endif

  CLI::App app{"Colvars stub interface for testing"};
  argv = app.ensure_utf8(argv);
  std::string configuration_file;
  std::string output_prefix;
  std::string trajectory_file;
  bool output_force = false;
  int trial_evaluations = 0;
  bool trial_consistency = false;
  app.add_option("-c,--configuration_file", configuration_file, "Input Colvars configuration file")->required(true);
  app.add_option("-t,--trajectory_file", trajectory_file, "Input trajectory file")->required(true);
  const auto output_prefix_option =
    app.add_option("-o,--output_prefix", output_prefix, "Output file prefix")->required(false);
  app.add_flag("--force,!--no-force", output_force, "Write the force files")->needs(output_prefix_option);
  app.add_option("--trial-evaluations", trial_evaluations,
                  "Call calc_trial() this many times on scaled coordinates before each calc()");
  app.add_flag("--trial-consistency", trial_consistency,
                "Require calc_trial() and calc() forces to match on the same coordinates");
  CLI11_PARSE(app, argc, argv);

  int err = 0;

  colvarproxy_stub *proxy = new colvarproxy_stub();
  colvarmodule *cvmodule = proxy->cvmodule;
  // Initialize simple unit system to test file input
  err |= proxy->set_unit_system("real", false);

  if (argc > 3) {
    err |= proxy->set_output_prefix(output_prefix);
  }
  err |= cvmodule->setup_input();
  err |= cvmodule->setup_output();
  err |= cvmodule->read_config_file(configuration_file.c_str());

  if (argc > 2) {
    // Read number of atoms from XYZ header
    std::ifstream ifs(trajectory_file);
    int natoms;
    ifs >> natoms;
    ifs.close();
    cvmodule->log("Reading trajectory for " + cvm::to_str(natoms)
              + " atoms from XYZ file " + trajectory_file);
    for (int ai = 0; ai < natoms; ai++) {
      proxy->init_atom(ai+1);
    }
    if ((trial_evaluations > 0) || trial_consistency) {
      std::string reasons;
      if (cvmodule->check_trial_support(reasons) != COLVARS_OK) {
        cvmodule->log("Trial evaluations are not supported by this configuration:\n" + reasons);
        // Tell compare_test.sh to skip this directory
        std::string const dir = output_prefix.substr(0, output_prefix.find_last_of('/'));
        std::ofstream(dir + "/disabled").close();
        delete proxy;
        return 77;
      }
    }
    int io_err = 0;
    while (!io_err) {
      io_err = proxy->read_frame_xyz(trajectory_file.c_str(), output_force,
                                    trial_evaluations, trial_consistency);
      if (cvmodule->get_error() != COLVARS_OK) { err |= cvmodule->get_error(); break; }
      if (!io_err) cvmodule->log("Frame " + cvm::to_str(cvmodule->step_absolute()));
    }
    proxy->post_run();
    cvmodule->log("Done");
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
}
