#include <iostream>
#include <fstream>

#include "colvargrid.h"
#include "colvargrid_integrate.h"
#include "colvarproxy.h"


int main (int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "\n\nOne argument needed: gradient multicol file name.\n";
    return 1;
  }

  colvarproxy *proxy = new colvarproxy();
  proxy->cvmodule = new colvarmodule(proxy); // This could be omitted if we used the colvarproxy_stub class

  std::string gradfile (argv[1]);
  std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile);
  if (proxy->cvmodule->get_error()) { return -1; }

  int itmax = 10000;
  cvm::real err;
  cvm::real tol = 1e-8;

  colvargrid_integrate fes(grad_ptr);
  fes.set_div();
  fes.integrate(itmax, tol, err);
  fes.set_zero_minimum();

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
