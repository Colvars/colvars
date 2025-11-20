#include <iostream>

#include "colvargrid.h"
#include "colvargrid_integrate.h"
#include "colvarproxy.h"

// Integrate provided gradients while monitoring convergence towards a provided scalar grid
// (typically the result of a previous integration)

int main (int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "\n\nOne argument needed: gradient multicol file name.\n";
    return 1;
  }

  colvarproxy *proxy = new colvarproxy();
  proxy->cvmodule = new colvarmodule(proxy); // This could be omitted if we used the colvarproxy_stub class

  std::string gradfile (argv[1]);
  std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile);
  if (cvmodule->get_error()) { return -1; }

  cvm::real err = 1.;
  cvm::real tol = 1e-10;

  colvargrid_integrate fes(grad_ptr);
  fes.set_div();

  // Load reference
  colvar_grid_scalar ref(gradfile + ".ref");
  if (cvmodule->get_error()) { return -1; }

  if (ref.number_of_points() != fes.number_of_points()) {
    cvmodule->error("Reference grid has wrong number of points: " + cvm::to_str(ref.number_of_points()) + "\n");
    return -1;
  }

  // Monitor convergence
  int rounds = 100;
  int steps_per_round = 10;

  std::cout << 0 << " " << 0 << " " << ref.grid_rmsd(fes) << std::endl;

  for (int i = 0; i < rounds && err > tol; i++) {
    fes.reset(0.);
    fes.integrate(steps_per_round * (i+1), tol, err);
    fes.set_zero_minimum();

    std::cout << (i+1)*steps_per_round << " " <<  err << " " << ref.grid_rmsd(fes) << std::endl;

    char buff[100];
    snprintf(buff, sizeof(buff), "%04i", steps_per_round * (i+1));
  }

  if (fes.num_variables() < 3) {
    std::cout << "\nWriting integrated free energy in multicol format to " + gradfile + ".int\n";
    fes.write_multicol(std::string(gradfile + ".int"), "integrated free energy");
  } else { // Write 3D grids to more convenient DX format
    std::cout << "\nWriting integrated free energy in OpenDX format to " + gradfile + ".int.dx\n";
    fes.write_opendx(std::string(gradfile + ".int.dx"), "integrated free energy");
  }

  delete proxy;
  return 0;
}
