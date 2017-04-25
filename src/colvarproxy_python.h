// -*- c++ -*-

#ifndef COLVARPROXY_PYTHON_H
#define COLVARPROXY_PYTHON_H

#include <vector>


/// Methods for using Python within Colvars
class colvarproxy_python {

public:

  /// Constructor
  colvarproxy_python();

  /// Destructor
  virtual ~colvarproxy_python();

  /// Is Python available? (trigger initialization if needed)
  int py_available();

  /// Python implementation of script_obj_to_str()
  char const *py_obj_to_str(unsigned char *obj);

  /// Get a PyObject pointer to the object
  void *py_get_object(void *parent_obj, char const *obj_name);

  /// Get a PyObject pointer to the function
  void *py_get_function(void *parent_obj, char const *func_name);

  /// Raise error if object does not exist
  void py_check_obj(void *parent_obj, char const *obj_name, void *py_obj);

  /// Wrap a STL vector into a NumPy one
  void *py_array_from_vector(std::vector<double> const &v, void *vector_obj);

  /// Create a new NumPy array
  void *py_new_array(size_t n);

  /// Get the C-array buffer of a NumPy array
  void *py_get_array_buffer(void *obj);

  /// Create and/or set a Python object representing a single colvarvalue
  void *py_set_colvarvalue(colvarvalue const &value,
                           void *value_obj);

  /// \brief Read a colvarvalue from a Python object
  int py_get_colvarvalue(void *value_obj, colvarvalue *value);

  /// \brief Create and/or set a Python list of colvarvalue objects
  void *py_set_colvarvalues(std::vector<colvarvalue const *> const &values,
                            void *values_obj);

  /// \brief Create and/or set a Python sequence of colvar gradients (Jacobian
  /// matrices); if all matrices are one-dimensional, it will try to create a
  /// NumPy vector for them, otherwise a list
  void *py_set_colvar_gradients(
          std::vector<cvm::matrix2d<cvm::real> > const &grads,
          void *grads_obj);

  /// \brief Read a set of colvar gradients from a Python object
  int py_get_colvar_gradients(
        void *grads_obj,
        std::vector<cvm::matrix2d<cvm::real> > *gradient);

  /// Calls the Python function calc_colvar_forces()
  int py_run_force_callback();

  /// \brief Calls the Python function func with arguments cvc_values,
  /// colvar_value (if not NULL) and *out_obj; if the result is not Py_None,
  /// it will be saved into *out_obj
  int py_run_colvar_callback(void *func, void *cvc_values, void *colvar_value,
                             void **out_obj);

protected:

  /// Pointer to Python main module
  void *py_main;
  /// Pointers to NumPy Module and functions
  void *py_np, *py_np_array, *py_np_zeros, *py_np_frombuffer;

  /// Set Python pointers (let main program call Py_Initialize())
  virtual void init_py_pointers();
};


#endif
