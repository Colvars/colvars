#include <iostream>

#include "colvarmodule.h"
#include "colvarvalue.h"


extern "C" int main(int argc, char *argv[]) {

  colvarvalue x(colvarvalue::type_unit3vector);
  x.rvector_value = cvm::rvector(1.0, 0.0, 0.0);
  x.apply_constraints();

  std::cout << "x               = " 
            << cvm::to_str(x, cvm::cv_width, cvm::cv_prec) << std::endl;

  {
    colvarvalue y(colvarvalue::type_unit3vector);
    y.rvector_value = cvm::rvector(0.0, 1.0, 0.0);
    y.apply_constraints();

    colvarvalue const xtoy = 0.5*x.dist2_grad(y);
    colvarvalue const ytox = 0.5*y.dist2_grad(x);
    std::cout << "y               = "
              << cvm::to_str(y, cvm::cv_width, cvm::cv_prec) << std::endl;
  
    std::cout << "x.dist2_grad(y) = " 
              << cvm::to_str(xtoy, cvm::cv_width, cvm::cv_prec) << std::endl;
    std::cout << "y.dist2_grad(x) = "
              << cvm::to_str(ytox, cvm::cv_width, cvm::cv_prec) << std::endl;
  }

  {
    colvarvalue y(colvarvalue::type_unit3vector);
    y.rvector_value = cvm::rvector(1.0, 1.0, 0.0);
    y.apply_constraints();

    colvarvalue const xtoy = 0.5*x.dist2_grad(y);
    colvarvalue const ytox = 0.5*y.dist2_grad(x);
    std::cout << "y               = "
              << cvm::to_str(y, cvm::cv_width, cvm::cv_prec) << std::endl;
  
    std::cout << "x.dist2_grad(y) = " 
              << cvm::to_str(xtoy, cvm::cv_width, cvm::cv_prec) << std::endl;
    std::cout << "y.dist2_grad(x) = "
              << cvm::to_str(ytox, cvm::cv_width, cvm::cv_prec) << std::endl;
  }

  // Below 10^-8, (1 - cos^2) = 1 in dist2_grad()
  for (size_t i = 0; i < 8; i++) {
    colvarvalue y(colvarvalue::type_unit3vector);
    cvm::real eps = cvm::integer_power(10.0, -1*i)/2.0;
    y.rvector_value = x.rvector_value + cvm::rvector(eps, 2.0*eps, -0.5*eps);
    y.apply_constraints();
    std::cout << std::endl;
    std::cout << "Epsilon    =    " << eps << std::endl;
    std::cout << "y          = "
              << cvm::to_str(y, cvm::cv_width, cvm::cv_prec) << std::endl;

    colvarvalue xtoy = 0.5*x.dist2_grad(y);

    colvarvalue xminusy = x - y;

    std::cout << "dist2_grad = " 
              << cvm::to_str(xtoy, cvm::cv_width, cvm::cv_prec) << std::endl;
    std::cout << "difference = " 
              << cvm::to_str(xminusy, cvm::cv_width, cvm::cv_prec) << std::endl;
    for (size_t c = 0; c < 3; c++) {
      xtoy[c] /= xminusy[c];
    }
    std::cout << "ratio      = " 
              << cvm::to_str(xtoy, cvm::cv_width, cvm::cv_prec) << std::endl;
  }

  return 0;
}
