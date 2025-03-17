#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include "colvargrid.h"
#include "colvarproxy.h"


int main (int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "One argument needed: file name.\n";
    return 1;
  }

  colvarproxy *proxy = new colvarproxy();
  proxy->colvars = new colvarmodule(proxy); // This may be omitted if we used the colvarproxy_stub class

  std::string gradfile (argv[1]);
  std::shared_ptr<colvar_grid_count> count_ptr;

  std::string countfile;

  // Look for matching count file
  size_t pos = gradfile.rfind(std::string(".czar.grad"));
  if (pos != std::string::npos) {
    countfile = gradfile.substr(0,pos) + ".zcount";
  } else {
    pos = gradfile.rfind(std::string(".grad"));
    if (pos != std::string::npos) {
      countfile = gradfile.substr(0,pos) + ".count";
    }
  }
  if (countfile.size()) {
    struct stat buffer;
    if (stat(countfile.c_str(), &buffer) == 0) {
      std::cout << "Found associated count file " << countfile << ", reading...\n";
      count_ptr.reset(new colvar_grid_count(countfile));
      if (!count_ptr || count_ptr->nd == 0) { // catch constructor failure
        cvm::error("Error reading count grid.");
        return cvm::get_error();
      }
    }
  }

  std::cout << "Reading gradient file " << gradfile << std::endl;
  std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile, count_ptr);
  if (!grad_ptr || grad_ptr->nd == 0) { // catch constructor failure
    cvm::error("Error reading gradient grid.");
    return cvm::get_error();
  }

  int itmax = 1000;
  cvm::real err;
  cvm::real tol = 1e-6;

  integrate_potential potential(grad_ptr);
  potential.set_div();
  potential.integrate(itmax, tol, err);
  potential.set_zero_minimum();

  std::cout << "Writing integrated potential file " << gradfile + ".int" << std::endl;
  potential.write_multicol(std::string(gradfile + ".int"), "integrated potential");

  std::cout << "Writing internal gradient to file " << gradfile + ".out" << std::endl;
  grad_ptr->write_multicol(std::string(gradfile + ".out"), "integrated potential");

  delete proxy;
  return 0;
}
