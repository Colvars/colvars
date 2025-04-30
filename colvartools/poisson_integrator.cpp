#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include "colvargrid.h"
#include "colvarproxy.h"


int main (int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "\n\nOne argument needed: gradient multicol file name.\n";
    return 1;
  }
  colvarproxy *proxy = new colvarproxy();
  proxy->colvars = new colvarmodule(proxy); // This could be omitted if we used the colvarproxy_stub class

  std::string gradfile (argv[1]);
  std::string countfile;
  std::shared_ptr<colvar_grid_count> count_ptr;

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

  count_ptr->write_multicol("counts.dat");

  std::cout << "Reading gradient file " << gradfile << std::endl;
  std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile, count_ptr);
  // std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile);

  if (!grad_ptr || grad_ptr->nd == 0) { // catch constructor failure
    cvm::error("Error reading gradient grid.");
    return cvm::get_error();
  }

  grad_ptr->write_multicol("gradient_in.dat");

  int itmax = 100;
  cvm::real err;
  cvm::real tol = 1e-8;

  integrate_potential potential(grad_ptr);
  potential.prepare_laplacian_calculation();
  // potential.print_laplacian_preparations();

  potential.print_laplacian_preparations();
  potential.set_div();
  // potential.set_weighted_div();

  colvar_grid_scalar div(potential);
  div.data = potential.divergence;
  div.write_multicol("divergence.dat");

  // std::vector<cvm::real> laplacian_matrix (potential.computation_grid->nt, 0);
  // std::vector<cvm::real> test_vector (potential.computation_grid->nt, 1);
  // std::vector<cvm::real> complete_div (potential.computation_grid->nt, 0);
  // //
  // potential.laplacian_weighted<true>(test_vector, laplacian_matrix);
  // for (int i = 0; i < potential.computation_grid->nt; i++){
  // complete_div[i] = potential.divergence[i] + potential.div_border_supplement[i];
  // }
  // std::cout << potential.divergence.size() << " " << potential.div_border_supplement.size() << std::endl;
  // saveVectorToCSV(complete_div, "divergence.csv");
  // // saveVectorToCSV(potential.laplacian_matrix_test, "laplacian.csv");


  potential.integrate(itmax, tol, err, true, true);
  potential.set_zero_minimum();

  if (potential.num_variables() < 3) {
    std::cout << "\nWriting integrated potential in multicol format to " + gradfile + ".int\n";
    potential.write_multicol(std::string(gradfile + ".int"), "integrated potential");
  } else { // Write 3D grids to more convenient DX format
    std::cout << "\nWriting integrated potential in OpenDX format to " + gradfile + ".int.dx\n";
    potential.write_opendx(std::string(gradfile + ".int.dx"), "integrated potential");
  }

  delete proxy;
  return 0;
}
