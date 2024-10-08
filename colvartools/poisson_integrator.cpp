#include <iostream>
#include <fstream>

#include "colvargrid.h"
#include "colvarproxy.h"


int main (int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "One argument needed: file name.\n";
    return 1;
  }

  colvarproxy *proxy = new colvarproxy();
  colvarmodule *colvars = new colvarmodule(proxy);

  std::string gradfile (argv[1]);
  std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile);

  int itmax = 1000;
  cvm::real err;
  cvm::real tol = 1e-6;

  integrate_potential potential(grad_ptr);
  potential.set_div();
  potential.integrate(itmax, tol, err);
  potential.set_zero_minimum();

  potential.write_multicol(std::string(gradfile + ".int"),
                           "integrated potential");

  delete colvars;
  return 0;
}
