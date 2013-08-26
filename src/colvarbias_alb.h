#ifndef COLVARBIAS_ALB_H
#define COLVARBIAS_ALB_H

#include "colvar.h"
#include "colvarbias.h"

class colvarbias_alb : public colvarbias {
  
public:
  colvarbias_alb(std::string const &conf, char const *key);
  virtual ~colvarbias_alb();
  

  virtual cvm::real update();

  /// Read the bias configuration from a restart file
  virtual std::istream & read_restart (std::istream &is);

  /// Write the bias configuration to a restart file
  virtual std::ostream & write_restart (std::ostream &os);

  /// Write a label to the trajectory file (comment line)
  virtual std::ostream & write_traj_label (std::ostream &os);

  /// Output quantities such as the bias energy to the trajectory file
  virtual std::ostream & write_traj (std::ostream &os);

protected:

  /// \brief Restraint centers
  std::vector<colvarvalue> colvar_centers;  

  /// \brief colvar moments, used for calculating the gradient
  std::vector<colvarvalue> means;
  std::vector<cvm::real> means_sq;
  std::vector<colvarvalue> means_cu;
  int update_calls;

  ///\brief how often to update coupling force
  int update_freq;

  ///\brief Estimated range of coupling force values
  cvm::real max_coupling_change;

  /// \brief accumated couping force; used in stochastic online gradient descent algorithm
  cvm::real coupling_force_accum;

  /// \brief current coupling force
  cvm::real coupling_force; 

  /// \brief equilibration time of the colvars
  int equil_time;

  // \brief if we're equilibrating our estimates or collecting data
  bool b_equilibration;

  /// \brief flag for outputting colvar centers
  bool b_output_centers;

  /// \brief flag for outputting current gradient
  bool b_output_grad;

  /// \brief flag for outputting coupling force
  bool b_output_coupling;

  cvm::real restraint_potential(cvm::real k,  colvar*  x, const colvarvalue& xcenter) const;

  /// \brief Force function
  colvarvalue restraint_force(cvm::real k,  colvar* x,  const colvarvalue& xcenter) const;

  ///\brief Unit scaling
  cvm::real restraint_convert_k(cvm::real k, cvm::real dist_measure) const;

};

#endif
