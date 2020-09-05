/// @file     jacobi_pd.h
/// @brief    Calculate the eigenvalues and eigevectors of a symmetric matrix
///           using the Jacobi eigenvalue algorithm.
/// @author   Andrew Jewett
/// @note     Benchmarks and tests: https://github.com/jewettaij/jacobi_pd
/// @note     A version of this code is included with LAMMPS ("math_eigen.h").
/// @license  CC0-1.0  (public domain)

#ifndef _MATH_EIGEN_H
#define _MATH_EIGEN_H

#include <algorithm>
#include <cmath>
//#include <cassert>


namespace MathEigen {

/// @brief  Allocate a 2-dimensional array.  (Uses row-major order.)
/// @note   (Matrix allocation functions can also be found in colvartypes.h)
template<typename Entry>
void Alloc2D(size_t nrows,          //!< size of the array (number of rows)
             size_t ncols,          //!< size of the array (number of columns)
             Entry ***paaX          //!< pointer to a 2D C-style array
             );

/// @brief  Deallocate arrays that were created using Alloc2D().
/// @note   (Matrix allocation functions can also be found in colvartypes.h)
template<typename Entry>
void Dealloc2D(Entry ***paaX        //!< pointer to a 2D C-style array
               );


// ---- Eigendecomposition of small dense symmetric matrices ----

/// @class Jacobi
/// @brief Calculate the eigenvalues and eigevectors of a symmetric matrix
///        using the Jacobi eigenvalue algorithm.  Code for the Jacobi class
///        (along with tests and benchmarks) is available free of copyright at
///        https://github.com/jewettaij/jacobi_pd
/// @note  The "Vector" and "Matrix" type arguments can be any
///        C or C++ object that support indexing, including pointers or vectors.
/// @details
///   Usage example:
/// @code
///
/// int n = 5;       // Matrix size
/// double **M;      // A symmetric n x n matrix you want to diagonalize
/// double *evals;   // Store the eigenvalues here.
/// double **evects; // Store the eigenvectors here.
/// // Allocate space for M, evals, and evects, and load contents of M (omitted)
///
/// // Now create an instance of Jacobi ("eigen_calc"). This will allocate space
/// // for storing intermediate calculations.  Once created, it can be reused
/// // multiple times without paying the cost of allocating memory on the heap.
///
/// Jacobi<double, double*, double**> eigen_calc(n);
///
/// // Note:
/// // If the matrix you plan to diagonalize (M) is read-only, use this instead:
/// // Jacobi<double, double*, double**, double const*const*> eigen_calc(n);
/// // If you prefer using vectors over C-style pointers, this works also:
/// // Jacobi<double, vector<double>&, vector<vector<double>>&> eigen_calc(n);
///
/// // Now, calculate the eigenvalues and eigenvectors of M
///
/// eigen_calc.Diagonalize(M, evals, evects);
///
/// @endcode

template<typename Scalar,
         typename Vector,
         typename Matrix,
         typename ConstMatrix=Matrix>
class Jacobi
{
  int n;            //!< the size of the matrix
  Scalar **M;       //!< local copy of the matrix being analyzed
  // Precomputed cosine, sine, and tangent of the most recent rotation angle:
  Scalar c;         //!< = cos(θ)
  Scalar s;         //!< = sin(θ)
  Scalar t;         //!< = tan(θ),  (note |t|<=1)
  int *max_idx_row; //!< for row i, the index j of the maximum element where j>i

public:

  /// @brief  Specify the size of the matrices you want to diagonalize later.
  /// @param n  the size (ie. number of rows) of the square matrix
  void SetSize(int n);

  Jacobi(int n = 0) {
    Init();
    SetSize(n);
  }

  ~Jacobi() {
    Dealloc();
  }

  // @typedef choose the criteria for sorting eigenvalues and eigenvectors
  typedef enum eSortCriteria {
    DO_NOT_SORT,
    SORT_DECREASING_EVALS,
    SORT_INCREASING_EVALS,
    SORT_DECREASING_ABS_EVALS,
    SORT_INCREASING_ABS_EVALS
  } SortCriteria;

  /// @brief Calculate the eigenvalues and eigevectors of a symmetric matrix
  ///        using the Jacobi eigenvalue algorithm.
  /// @returns 0 if the algorithm converged,
  ///          1 if the algorithm failed to converge. (IE, if the number of
  ///            pivot iterations exceeded max_num_sweeps * iters_per_sweep,
  ///            where iters_per_sweep = (n*(n-1)/2))
  /// @note  To reduce the computation time further, set calc_evecs=false.
  int
  Diagonalize(ConstMatrix mat, //!< the matrix you wish to diagonalize (size n)
              Vector eval,   //!< store the eigenvalues here
              Matrix evec,   //!< store the eigenvectors here (in rows)
              SortCriteria sort_criteria=SORT_DECREASING_EVALS,//!<sort results?
              bool calc_evecs=true,      //!< calculate the eigenvectors?
              int max_num_sweeps = 50);  //!< limit the number of iterations

private:

  /// @brief Calculate the components of a rotation matrix which performs a
  ///        rotation in the i,j plane by an angle (θ) causing M[i][j]=0.
  ///        The resulting parameters will be stored in this->c, this->s, and
  ///        this->t (which store cos(θ), sin(θ), and tan(θ), respectively).
  void CalcRot(Scalar const *const *M,   //!< matrix
               int i,      //!< row index
               int j);     //!< column index

  /// @brief Apply the (previously calculated) rotation matrix to matrix M
  ///        by multiplying it on both sides (a "similarity transform").
  ///        (To save time, only the elements in the upper-right-triangular
  ///         region of the matrix are updated.  It is assumed that i < j.)
  void ApplyRot(Scalar **M,  //!< matrix
                int i,     //!< row index
                int j);    //!< column index

  /// @brief Multiply matrix E on the left by the (previously calculated)
  ///        rotation matrix.
  void ApplyRotLeft(Matrix E,  //!< matrix
                    int i,     //!< row index
                    int j);    //!< column index

  ///@brief Find the off-diagonal index in row i whose absolute value is largest
  int MaxEntryRow(Scalar const *const *M, int i) const;

  /// @brief Find the indices (i_max, j_max) marking the location of the
  ///        entry in the matrix with the largest absolute value.  This
  ///        uses the max_idx_row[] array to find the answer in O(n) time.
  /// @returns This function does not return a avalue.  However after it is
  ///          invoked, the location of the largest matrix element will be
  ///          stored in the i_max and j_max arguments.
  void MaxEntry(Scalar const *const *M, int& i_max, int& j_max) const;

  // @brief Sort the rows in M according to the numbers in v (also sorted)
  void SortRows(Vector v, //!< vector containing the keys used for sorting
                Matrix M, //!< matrix whose rows will be sorted according to v
                int n,    //!< size of the vector and matrix
                SortCriteria s=SORT_DECREASING_EVALS //!< sort decreasing order?
                ) const;

  // memory management:
  void Alloc(int n);
  void Init();
  void Dealloc();

public:
  // memory management: copy and move constructor, swap, and assignment operator
  Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMatrix>& source);
  Jacobi(Jacobi<Scalar, Vector, Matrix, ConstMatrix>&& other);
  void swap(Jacobi<Scalar, Vector, Matrix, ConstMatrix> &other);
  Jacobi<Scalar, Vector, Matrix, ConstMatrix>& operator = (Jacobi<Scalar, Vector, Matrix, ConstMatrix> source);

}; // class Jacobi





// -------------- IMPLEMENTATION --------------

// Utilities for memory allocation

template<typename Entry>
void Alloc2D(size_t nrows,          // size of the array (number of rows)
             size_t ncols,          // size of the array (number of columns)
             Entry ***paaX)         // pointer to a 2D C-style array
{
  //assert(paaX);
  *paaX = new Entry* [nrows];  //conventional 2D C array (pointer-to-pointer)
  (*paaX)[0] = new Entry [nrows * ncols];  // 1D C array (contiguous memory)
  for(size_t iy=0; iy<nrows; iy++)
    (*paaX)[iy] = (*paaX)[0] + iy*ncols;
  // The caller can access the contents using (*paaX)[i][j]
}

template<typename Entry>
void Dealloc2D(Entry ***paaX)       // pointer to a 2D C-style array
{
  if (paaX && *paaX) {
    delete [] (*paaX)[0];
    delete [] (*paaX);
    *paaX = nullptr;
  }
}



template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Diagonalize(ConstMatrix mat,    // the matrix you wish to diagonalize (size n)
            Vector eval,        // store the eigenvalues here
            Matrix evec,        // store the eigenvectors here (in rows)
            SortCriteria sort_criteria, // sort results?
            bool calc_evec,     // calculate the eigenvectors?
            int max_num_sweeps) // limit the number of iterations ("sweeps")
{
  // -- Initialization --
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)          //copy mat[][] into M[][]
      M[i][j] = mat[i][j];               //(M[][] is a local copy we can modify)

  if (calc_evec)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        evec[i][j] = (i==j) ? 1.0 : 0.0; //Set evec equal to the identity matrix

  for (int i = 0; i < n-1; i++)          //Initialize the "max_idx_row[]" array 
    max_idx_row[i] = MaxEntryRow(M, i);  //(which is needed by MaxEntry())

  // -- Iteration --
  int n_iters;
  int max_num_iters = max_num_sweeps*n*(n-1)/2; //"sweep" = n*(n-1)/2 iters
  for (n_iters=0; n_iters < max_num_iters; n_iters++) {
    int i,j;
    MaxEntry(M, i, j); // Find the maximum entry in the matrix. Store in i,j

    // If M[i][j] is small compared to M[i][i] and M[j][j], set it to 0.
    if ((M[i][i] + M[i][j] == M[i][i]) && (M[j][j] + M[i][j] == M[j][j])) {
      M[i][j] = 0.0;
      max_idx_row[i] = MaxEntryRow(M,i); //must also update max_idx_row[i]
    }

    if (M[i][j] == 0.0)
      break;

    // Otherwise, apply a rotation to make M[i][j] = 0
    CalcRot(M, i, j);  // Calculate the parameters of the rotation matrix.
    ApplyRot(M, i, j); // Apply this rotation to the M matrix.
    if (calc_evec)     // Optional: If the caller wants the eigenvectors, then
      ApplyRotLeft(evec,i,j); // apply the rotation to the eigenvector matrix

  } //for (int n_iters=0; n_iters < max_num_iters; n_iters++)

  // -- Post-processing --
  for (int i = 0; i < n; i++)
    eval[i] = M[i][i];

  // Optional: Sort results by eigenvalue.
  SortRows(eval, evec, n, sort_criteria);

  return (n_iters == max_num_iters);
}


/// @brief Calculate the components of a rotation matrix which performs a
///        rotation in the i,j plane by an angle (θ) that (when multiplied on
///        both sides) will zero the ij'th element of M, so that afterwards
///        M[i][j] = 0.  The results will be stored in c, s, and t
///        (which store cos(θ), sin(θ), and tan(θ), respectively).

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
CalcRot(Scalar const *const *M,    // matrix
        int i,       // row index
        int j)       // column index
{
  t = 1.0; // = tan(θ)
  Scalar M_jj_ii = (M[j][j] - M[i][i]);
  if (M_jj_ii != 0.0) {
    // kappa = (M[j][j] - M[i][i]) / (2*M[i][j])
    Scalar kappa = M_jj_ii;
    t = 0.0;
    Scalar M_ij = M[i][j];
    if (M_ij != 0.0) {
      kappa /= (2.0*M_ij);
      // t satisfies: t^2 + 2*t*kappa - 1 = 0
      // (choose the root which has the smaller absolute value)
      t = 1.0 / (std::sqrt(1 + kappa*kappa) + std::abs(kappa));
      if (kappa < 0.0)
        t = -t;
    }
  }
  //assert(std::abs(t) <= 1.0);
  c = 1.0 / std::sqrt(1 + t*t);
  s = c*t;
}

/// brief  Perform a similarity transformation by multiplying matrix M on both
///         sides by a rotation matrix (and its transpose) to eliminate M[i][j].
/// details This rotation matrix performs a rotation in the i,j plane by
///         angle θ.  This function assumes that c=cos(θ). s=som(θ), t=tan(θ)
///         have been calculated in advance (using the CalcRot() function).
///         It also assumes that i<j.  The max_idx_row[] array is also updated.
///         To save time, since the matrix is symmetric, the elements
///         below the diagonal (ie. M[u][v] where u>v) are not computed.
/// verbatim
///   M' = R^T * M * R
/// where R the rotation in the i,j plane and ^T denotes the transpose.
///                 i         j
///       _                             _
///      |  1                            | 
///      |    .                          |
///      |      .                        |
///      |        1                      |
///      |          c   ...   s          |
///      |          .  .      .          |
/// R  = |          .    1    .          |
///      |          .      .  .          |
///      |          -s  ...   c          |
///      |                      1        |
///      |                        .      |
///      |                          .    |
///      |_                           1 _|
/// endverbatim
///
/// Let M' denote the matrix M after multiplication by R^T and R.
/// The components of M' are:
///
/// verbatim
///   M'_uv =  Σ_w  Σ_z   R_wu * M_wz * R_zv
/// endverbatim
///
/// Note that a the rotation at location i,j will modify all of the matrix
/// elements containing at least one index which is either i or j
/// such as: M[w][i], M[i][w], M[w][j], M[j][w].
/// Check and see whether these modified matrix elements exceed the 
/// corresponding values in max_idx_row[] array for that row.
/// If so, then update max_idx_row for that row.
/// This is somewhat complicated by the fact that we must only consider
/// matrix elements in the upper-right triangle strictly above the diagonal.
/// (ie. matrix elements whose second index is > the first index).
/// The modified elements we must consider are marked with an "X" below:
///
/// @verbatim
///                 i         j
///       _                             _
///      |  .       X         X          | 
///      |    .     X         X          |
///      |      .   X         X          |
///      |        . X         X          |
///      |          X X X X X 0 X X X X  |  i
///      |            .       X          |
///      |              .     X          |
/// M  = |                .   X          |
///      |                  . X          |
///      |                    X X X X X  |  j
///      |                      .        |
///      |                        .      |
///      |                          .    |
///      |_                           . _|
/// @endverbatim

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
ApplyRot(Scalar **M,  // matrix
         int i,     // row index
         int j)     // column index
{
  // Recall that:
  // c = cos(θ)
  // s = sin(θ)
  // t = tan(θ) (which should be <= 1.0)

  // Compute the diagonal elements of M which have changed:
  M[i][i] -= t * M[i][j];
  M[j][j] += t * M[i][j];
  // Note: This is algebraically equivalent to:
  // M[i][i] = c*c*M[i][i] + s*s*M[j][j] - 2*s*c*M[i][j]
  // M[j][j] = s*s*M[i][i] + c*c*M[j][j] + 2*s*c*M[i][j]

  //Update the off-diagonal elements of M which will change (above the diagonal)

  //assert(i < j);

  M[i][j] = 0.0;

  //compute M[w][i] and M[i][w] for all w!=i,considering above-diagonal elements
  for (int w=0; w < i; w++) {        // 0 <= w <  i  <  j < n
    M[i][w] = M[w][i]; //backup the previous value. store below diagonal (i>w)
    M[w][i] = c*M[w][i] - s*M[w][j]; //M[w][i], M[w][j] from previous iteration
    if (i == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][i])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=i;
    //assert(max_idx_row[w] == MaxEntryRow(M, w));
  }
  for (int w=i+1; w < j; w++) {      // 0 <= i <  w  <  j < n
    M[w][i] = M[i][w]; //backup the previous value. store below diagonal (w>i)
    M[i][w] = c*M[i][w] - s*M[w][j]; //M[i][w], M[w][j] from previous iteration
  }
  for (int w=j+1; w < n; w++) {      // 0 <= i < j+1 <= w < n
    M[w][i] = M[i][w]; //backup the previous value. store below diagonal (w>i)
    M[i][w] = c*M[i][w] - s*M[j][w]; //M[i][w], M[j][w] from previous iteration
  }

  // now that we're done modifying row i, we can update max_idx_row[i]
  max_idx_row[i] = MaxEntryRow(M, i);

  //compute M[w][j] and M[j][w] for all w!=j,considering above-diagonal elements
  for (int w=0; w < i; w++) {        // 0 <=  w  <  i <  j < n
    M[w][j] = s*M[i][w] + c*M[w][j]; //M[i][w], M[w][j] from previous iteration
    if (j == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][j])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=j;
    //assert(max_idx_row[w] == MaxEntryRow(M, w));
  }
  for (int w=i+1; w < j; w++) {      // 0 <= i+1 <= w <  j < n
    M[w][j] = s*M[w][i] + c*M[w][j]; //M[w][i], M[w][j] from previous iteration
    if (j == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][j])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=j;
    //assert(max_idx_row[w] == MaxEntryRow(M, w));
  }
  for (int w=j+1; w < n; w++) {      // 0 <=  i  <  j <  w < n
    M[j][w] = s*M[w][i] + c*M[j][w]; //M[w][i], M[j][w] from previous iteration
  }
  // now that we're done modifying row j, we can update max_idx_row[j]
  max_idx_row[j] = MaxEntryRow(M, j);

} //Jacobi::ApplyRot()




///@brief Multiply matrix M on the LEFT side by a transposed rotation matrix R^T
///       This matrix performs a rotation in the i,j plane by angle θ  (where
///       the arguments "s" and "c" refer to cos(θ) and sin(θ), respectively).
/// @verbatim
///   E'_uv = Σ_w  R_wu * E_wv
/// @endverbatim

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
ApplyRotLeft(Matrix E,  // matrix
             int i,     // row index
             int j)     // column index
{
  // Recall that c = cos(θ) and s = sin(θ)
  for (int v = 0; v < n; v++) {
    Scalar Eiv = E[i][v]; //backup E[i][v]
    E[i][v] = c*E[i][v] - s*E[j][v];
    E[j][v] = s*Eiv     + c*E[j][v];
  }
}



template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
MaxEntryRow(Scalar const *const *M, int i) const {
  int j_max = i+1;
  for(int j = i+2; j < n; j++)
    if (std::abs(M[i][j]) > std::abs(M[i][j_max]))
      j_max = j;
  return j_max;
}



template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
MaxEntry(Scalar const *const *M, int& i_max, int& j_max) const {
  // find the maximum entry in the matrix M in O(n) time
  i_max = 0;
  j_max = max_idx_row[i_max];
  Scalar max_entry = std::abs(M[i_max][j_max]);
  int nm1 = n-1;
  for (int i=1; i < nm1; i++) {
    int j = max_idx_row[i];
    if (std::abs(M[i][j]) > max_entry) {
      max_entry = std::abs(M[i][j]);
      i_max = i;
      j_max = j;
    }
  }
  //#ifndef NDEBUG
  //// make sure that the maximum element really is stored at i_max, j_max
  ////   (comment out the next loop before using. it slows down the code)
  //for (int i = 0; i < nm1; i++)
  //  for (int j = i+1; j < n; j++)
  //    assert(std::abs(M[i][j]) <= max_entry);
  //#endif
}



//Sort the rows of a matrix "evec" by the numbers contained in "eval"
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
SortRows(Vector eval, Matrix evec, int n, SortCriteria sort_criteria) const
{
  for (int i = 0; i < n-1; i++) {
    int i_max = i;
    for (int j = i+1; j < n; j++) {
      // find the "maximum" element in the array starting at position i+1
      switch (sort_criteria) {
      case SORT_DECREASING_EVALS:
        if (eval[j] > eval[i_max])
          i_max = j;
        break;
      case SORT_INCREASING_EVALS:
        if (eval[j] < eval[i_max])
          i_max = j;
        break;
      case SORT_DECREASING_ABS_EVALS:
        if (std::abs(eval[j]) > std::abs(eval[i_max]))
          i_max = j;
        break;
      case SORT_INCREASING_ABS_EVALS:
        if (std::abs(eval[j]) < std::abs(eval[i_max]))
          i_max = j;
        break;
      default:
        break;
      }
    }
    std::swap(eval[i], eval[i_max]); // sort "eval"
    for (int k = 0; k < n; k++)
      std::swap(evec[i][k], evec[i_max][k]); // sort "evec"
  }
}



template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Init() {
  n = 0;
  M = nullptr;
  max_idx_row = nullptr;
}

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
SetSize(int n) {
  Dealloc();
  Alloc(n);
}

// memory management:

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Alloc(int n) {
  this->n = n;
  if (n > 0) {
    max_idx_row = new int[n];
    Alloc2D(n, n, &M);
  }
}

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Dealloc() {
  if (max_idx_row)
    delete [] max_idx_row;
  Dealloc2D(&M);
  Init();
}

// memory management: copy and move constructor, swap, and assignment operator:

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMatrix>& source)
{
  Init();
  SetSize(source.n);
  //assert(n == source.n);

  // The following lines aren't really necessary, because the contents
  // of source.M and source.max_idx_row are not needed (since they are
  // overwritten every time Jacobi::Diagonalize() is invoked).
  std::copy(source.max_idx_row,
            source.max_idx_row + n,
            max_idx_row);
  for (int i = 0; i < n; i++)
    std::copy(source.M[i],
              source.M[i] + n,
              M[i]);
}

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
swap(Jacobi<Scalar, Vector, Matrix, ConstMatrix> &other) {
  std::swap(n, other.n);
  std::swap(max_idx_row, other.max_idx_row);
  std::swap(M, other.M);
}

// Move constructor (C++11)
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Jacobi(Jacobi<Scalar, Vector, Matrix, ConstMatrix>&& other) {
  Init();
  this->swap(other);
}

// Using the "copy-swap" idiom for the assignment operator
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>&
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
operator = (Jacobi<Scalar, Vector, Matrix, ConstMatrix> source) {
  this->swap(source);
  return *this;
}

} // namespace MathEigen


#endif //#ifndef _MATH_EIGEN_H

