#include <fstream>

#include "colvargrid.h"
#include "colvarproxy.h"


int main (int argc, char *argv[]) {

  colvarproxy *proxy = new colvarproxy();
  colvarmodule *colvars = new colvarmodule(proxy);

  if (argc < 2) {
    std::cerr << "One argument needed: file name.\n";
    return 1;
  }

  std::string gradfile (argv[1]);
  colvar_grid_gradient grad(gradfile);

  int itmax = 1000;
  cvm::real err;
  cvm::real tol = 1e-6;

  integrate_potential potential(&grad);
  potential.set_div();
  potential.integrate(itmax, tol, err);
  potential.set_zero_minimum();

  std::cout << "Writing integrated potential to " + gradfile + ".int\n";
  std::ofstream os(std::string(gradfile + ".int").c_str());
  potential.write_multicol(os);
  os.close();

  return 0;
}