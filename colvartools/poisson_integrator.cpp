#include <iostream>
#include <fstream>

#include "colvargrid.h"
#include "colvarproxy.h"


int main (int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "\n\nOne argument needed: gradient multicol file name.\n";
    return 1;
  }

  colvarproxy *proxy = new colvarproxy();
  colvarmodule *colvars = new colvarmodule(proxy);

  std::string gradfile (argv[1]);
  std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile);
  if (cvm::get_error()) { return -1; }

  int itmax = 10000;
  cvm::real err;
  cvm::real tol = 1e-8;

  integrate_potential potential(grad_ptr);
  potential.set_div();
  potential.integrate(itmax, tol, err);
  potential.set_zero_minimum();

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
