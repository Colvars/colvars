// -*- c++ -*-

#ifndef COLVARCOMP_DEER_H
#define COLVARCOMP_DEER_H


class colvar;

/// \brief DEER time-signal colvar component
/// (Marinelli and Fiorin, Structure, 2018)
class colvar::deer_kernel
  : public colvar::cvc
{
protected:

  /// kernel calculation routine
  double kdeer(cvm::real const &r, cvm::real const &t);
  /// deer time trace file
  std::string deer_time_file;
  /// size of deer file
  int deersize;
  /// times vector
  cvm::vector1d<cvm::real> timesdeer;
  /// experimental deer values
  cvm::vector1d<cvm::real> deerexpvalues;
  /// width for deer grid
  cvm::real                deerwidth;
  /// lower boundary for deer grid
  cvm::real                deerlower;
  /// upper boundary for deer grid
  cvm::real                deerupper;
  /// number of points in discretized deer kernel function
  int rpoints;
  /// discretized deer kernel function
  cvm::matrix2d<cvm::real> deerk;
  /// First atom group
  cvm::atom_group  *group1;
  /// Second atom group
  cvm::atom_group  *group2;
  /// Vector distance, cached to be recycled
  cvm::rvector     dist_v;
  /// Derivatives vector , cached to be recycled
  cvm::vector1d<cvm::real> deerder;

public:

  deer_kernel(std::string const &conf);
  deer_kernel();
  virtual ~deer_kernel() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;
};


#endif
