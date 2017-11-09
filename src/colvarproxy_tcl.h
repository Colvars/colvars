// -*- c++ -*-

#ifndef COLVARPROXY_TCL_H
#define COLVARPROXY_TCL_H

#include <vector>


/// Methods for using Tcl within Colvars
class colvarproxy_tcl {

public:

  /// Constructor
  colvarproxy_tcl();

  /// Destructor
  virtual ~colvarproxy_tcl();

  /// Is Tcl available? (trigger initialization if needed)
  int tcl_available();

  /// Tcl implementation of script_obj_to_str()
  char const *tcl_get_str(void *obj);

  // TODO: support for 64-bit integers is deferred until a consistent protocol
  // for treating them is in place in VMD and NAMD

  /// Get an integer from a Tcl object
  int tcl_get_int(void *obj, int *result);

  /// Get a real number from a Tcl object
  int tcl_get_real(void *obj, cvm::real *result);

  /// Get a vector of integers from a Tcl object
  int tcl_get_int_vector(void *obj, std::vector<int> *vec);

  /// Get a vector of reals from a Tcl object
  int tcl_get_real_vector(void *obj, std::vector<cvm::real> *vec);

  /// Copy a string into a Tcl object
  void *tcl_set_string(std::string const &s, void *obj);

  /// Assign a real number to a Tcl object
  void *tcl_set_int(int const &x, void *obj);

  /// Assign a real number to a Tcl object
  void *tcl_set_real(cvm::real const &x, void *obj);

  /// Assign a vector of integers to a Tcl list
  void *tcl_list_from_int_vec(std::vector<int> const &x);

  /// Assign a vector of reals to a Tcl list
  void *tcl_list_from_real_vec(std::vector<cvm::real> const &x);

  /// Assign a vector of 3-dimensional vectors to a Tcl list
  void *tcl_list_from_rvector_vec(std::vector<cvm::rvector> const &x);

  /// Tcl implementation of run_force_callback()
  int tcl_run_force_callback();

  /// Tcl implementation of run_colvar_callback()
  int tcl_run_colvar_callback(
              std::string const &name,
              std::vector<const colvarvalue *> const &cvcs,
              colvarvalue &value);

  /// Tcl implementation of run_colvar_gradient_callback()
  int tcl_run_colvar_gradient_callback(
              std::string const &name,
              std::vector<const colvarvalue *> const &cvcs,
              std::vector<cvm::matrix2d<cvm::real> > &gradient);

  /// Get a pointer to the Tcl interpreter
  inline void *get_tcl_interp()
  {
    return tcl_interp_;
  }

protected:

  /// Pointer to Tcl interpreter object
  void *tcl_interp_;

  /// Set Tcl pointers
  virtual void init_tcl_pointers();
};


#endif
