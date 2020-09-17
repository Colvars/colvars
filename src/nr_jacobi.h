// -*- c++ -*-

// TODO These functions will be removed in the near future and replaced wtih
// LGPL-compatible functions

namespace NR_Jacobi {

  /// Numerical recipes diagonalization
  int jacobi(cvm::real **a, cvm::real *d, cvm::real **v, int *nrot);

  /// Eigenvector sort
  int eigsrt(cvm::real *d, cvm::real **v);

  /// Transpose the matrix
  int transpose(cvm::real **v);

}

