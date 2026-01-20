// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.
//

#ifndef COLVARCOMP_COORDNUM_H
#define COLVARCOMP_COORDNUM_H

#include <memory>

#include "colvar.h"
#include "colvarcomp.h"
#include "colvarmodule.h"


/// \brief Colvar component: coordination number between two groups
/// (colvarvalue::type_scalar type, range [0:N1*N2])
class colvar::coordnum : public colvar::cvc {
public:

  coordnum();
  virtual ~coordnum();
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();

  enum {
    ef_null = 0,
    ef_gradients = 1,
    ef_anisotropic = (1 << 8),
    ef_use_pairlist = (1 << 9),
    ef_rebuild_pairlist = (1 << 10)
  };

  /// \brief Calculate a coordination number through the function
  /// (1-x**n)/(1-x**m), where x = |A1-A2|/r0 \param r0, r0_vec "cutoff" for
  /// the coordination number (scalar or vector depending on user choice)
  /// \param en Numerator exponent \param ed Denominator exponent \param First
  /// atom \param Second atom \param pairlist_elem pointer to pair flag for
  /// this pair \param tolerance A pair is defined as having a larger
  /// coordination than this number
  template <int flags>
  static cvm::real switching_function(cvm::real const &r0, cvm::rvector const &inv_r0_vec,
                                      cvm::rvector const &inv_r0sq_vec, int en, int ed,
                                      const cvm::real a1x, const cvm::real a1y, const cvm::real a1z,
                                      const cvm::real a2x, const cvm::real a2y, const cvm::real a2z,
                                      cvm::real &g1x, cvm::real &g1y, cvm::real &g1z,
                                      cvm::real &g2x, cvm::real &g2y, cvm::real &g2z,
                                      bool **pairlist_elem, cvm::real tolerance);

  /// Workhorse function
  template <int flags> int compute_coordnum();

  /// Workhorse function
  template <int flags> void main_loop(bool **pairlist_elem);

protected:
  /// First atom group
  cvm::atom_group *group1 = nullptr;
  /// Second atom group
  cvm::atom_group *group2 = nullptr;
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real r0;
  /// \brief "Cutoff vector" for anisotropic calculation
  cvm::rvector r0_vec;
  /// \brief Whether r/r0 or \vec{r}*\vec{1/r0_vec} should be used
  bool b_anisotropic = false;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 12;

  /// If true, group2 will be treated as a single atom
  bool b_group2_center_only = false;

  /// Tolerance for the pair list
  cvm::real tolerance = 0.0;

  /// Frequency of update of the pair list
  int pairlist_freq = 100;

  /// Pair list
  bool *pairlist = nullptr;

};


/// \brief Colvar component: self-coordination number within a group
/// (colvarvalue::type_scalar type, range [0:N*(N-1)/2])
class colvar::selfcoordnum : public colvar::cvc {
public:

  selfcoordnum();
  ~selfcoordnum();
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();

  /// Main workhorse function
  template <int flags> int compute_selfcoordnum();

protected:
  /// Selected atoms
  cvm::atom_group *group1 = nullptr;
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real r0;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 12;
  cvm::real tolerance = 0.0;
  int pairlist_freq = 100;

  bool *pairlist = nullptr;

};


/// \brief Colvar component: coordination number between two groups
/// (colvarvalue::type_scalar type, range [0:N1*N2])
class colvar::groupcoordnum : public colvar::distance {
public:
  groupcoordnum();
  virtual ~groupcoordnum() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();

protected:
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real r0;
  /// \brief "Cutoff vector" for anisotropic calculation
  cvm::rvector r0_vec;
  /// \brief Wheter dist/r0 or \vec{dist}*\vec{1/r0_vec} should ne be
  /// used
  bool b_anisotropic = false;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 12;
};


/// \brief Colvar component: hydrogen bond, defined as the product of
/// a colvar::coordnum and 1/2*(1-cos((180-ang)/ang_tol))
/// (colvarvalue::type_scalar type, range [0:1])
class colvar::h_bond : public colvar::cvc {
public:
  /// Constructor for atoms already allocated
  h_bond(cvm::atom_group::simple_atom const &acceptor, cvm::atom_group::simple_atom const &donor,
         cvm::real r0, int en, int ed);
  h_bond();
  virtual ~h_bond() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();

protected:
  /// \brief "Cutoff" distance between acceptor and donor
  cvm::real r0;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 8;
};


#endif // COLVARCOMP_COORDNUM_H
