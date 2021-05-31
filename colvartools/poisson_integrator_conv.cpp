#include <iostream>

#include "colvargrid.h"
#include "colvarproxy.h"

// Integrate provided gradients while monitoring convergence towards a provided scalar grid
// (typically the result of a previous integration)

int main (int argc, char *argv[])
{
  colvarproxy *proxy = new colvarproxy();
  colvarmodule *colvars = new colvarmodule(proxy);

  if (argc < 2) {
    std::cerr << "\n\nOne argument needed: gradient multicol file name.\n";
    return 1;
  }

  std::string gradfile (argv[1]);
  std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile);
  if (cvm::get_error()) { return -1; }

  cvm::real err = 1.;
  cvm::real tol = 1e-10;

  integrate_potential potential(grad_ptr);
  potential.set_div();

  // Load reference
  colvar_grid_scalar ref(gradfile + ".ref");
  if (cvm::get_error()) { return -1; }

  if (ref.number_of_points() != potential.number_of_points()) {
    cvm::error("Reference grid has wrong number of points: " + cvm::to_str(ref.number_of_points()) + "\n");
    return -1;
  }

  // Monitor convergence
  int rounds = 100;
  int steps_per_round = 10;

  std::cout << 0 << " " << 0 << " " << ref.grid_rmsd(potential) << std::endl;

  for (int i = 0; i < rounds && err > tol; i++) {
    potential.reset(0.);
    potential.integrate(steps_per_round * (i+1), tol, err);
    potential.set_zero_minimum();

    std::cout << (i+1)*steps_per_round << " " <<  err << " " << ref.grid_rmsd(potential) << std::endl;

    char buff[100];
    snprintf(buff, sizeof(buff), "%04i", steps_per_round * (i+1));
  }

  if (potential.num_variables() < 3) {
    std::cout << "\nWriting integrated potential in multicol format to " + gradfile + ".int\n";
    potential.write_multicol(std::string(gradfile + ".int"), "integrated potential");
  } else { // Write 3D grids to more convenient DX format
    std::cout << "\nWriting integrated potential in OpenDX format to " + gradfile + ".int.dx\n";
    potential.write_opendx(std::string(gradfile + ".int.dx"), "integrated potential");
  }
  delete colvars;
  return 0;
}
