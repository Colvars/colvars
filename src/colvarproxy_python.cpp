// -*- c++ -*-

#if defined(LMP_PYTHON) || defined(NAMD_PYTHON) || defined(VMDPYTHON)
#define COLVARS_PYTHON
#include "Python.h"
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_python.h"



colvarproxy_python::colvarproxy_python()
{
  py_main = NULL;
  py_np = py_np_array = py_np_zeros = py_np_frombuffer = NULL;
}


colvarproxy_python::~colvarproxy_python()
{
#if defined(COLVARS_PYTHON)
  Py_XDECREF(py_main);
  Py_XDECREF(py_np);
  Py_XDECREF(py_np_array);
  Py_XDECREF(py_np_zeros);
  Py_XDECREF(py_np_frombuffer);
#endif
}


int colvarproxy_python::py_available()
{
#if defined(COLVARS_PYTHON)
  if (py_main == NULL) init_py_pointers();
  if (cvm::debug()) {
    cvm::log("Checking if Python is available, returning yes.\n");
  }
  return COLVARS_OK;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


void colvarproxy_python::init_py_pointers()
{
#if defined(COLVARS_PYTHON)
  if (!Py_IsInitialized()) {
    cvm::error("Error: Python is not initialized.\n", BUG_ERROR);
    return;
  }

  PyObject *py_modules = PyImport_GetModuleDict(); // borrowed
  if (py_modules == NULL) {
    cvm::error("Error: cannot get list of Python modules.\n",
               BUG_ERROR);
    return;
  }

  if (PyMapping_HasKeyString(py_modules, const_cast<char *>("__main__"))) {
    py_main = PyMapping_GetItemString(py_modules,
                                      const_cast<char *>("__main__"));
  } else {
    cvm::error("Error: cannot get Python module __main__.\n", BUG_ERROR);
    return;
  }

  if (PyMapping_HasKeyString(py_modules, const_cast<char *>("numpy"))) {
    py_np = py_get_object(NULL, "numpy");
    py_np_array = py_get_object(py_np, "array");
    py_np_zeros = py_get_object(py_np, "zeros");
    py_np_frombuffer = py_get_object(py_np, "frombuffer");
  }
#endif
}


char const *colvarproxy_python::py_obj_to_str(unsigned char *obj)
{
#if defined(COLVARS_PYTHON)
  PyObject *po = reinterpret_cast<PyObject *>(obj);
  PyObject *po_repr = PyObject_Repr(po);
#if (PY_MAJOR_VERSION == 2)
  char const *str = PyString_AsString(po_repr);
#elif (PY_MAJOR_VERSION == 3)
  char const *str =  PyUnicode_AsUTF8(po);
#else
  "Unsupported Python?";
#endif
  Py_XDECREF(po_repr);
  return str;
#else
  return NULL;
#endif
}


void *colvarproxy_python::py_get_object(void *parent_obj,
                                        char const *obj_name)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_parent_obj =
    reinterpret_cast<PyObject *>(parent_obj ? parent_obj : py_main);
  if (py_main == NULL) init_py_pointers();
  PyGILState_STATE py_gil_state = PyGILState_Ensure();
  PyObject *py_obj = PyObject_GetAttrString(py_parent_obj, obj_name);
  PyGILState_Release(py_gil_state);
  return py_obj;
#else
  return NULL;
#endif
}


void *colvarproxy_python::py_get_function(void *parent_obj,
                                          char const *func_name)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_func =
    reinterpret_cast<PyObject *>(py_get_object(parent_obj, func_name));
  if (py_func != NULL) {
    if (!PyCallable_Check(py_func)) {
      cvm::error("Error: Python attribute "+
                 std::string(func_name)+" is not callable.\n",
                 INPUT_ERROR);
    }
  }
  return py_func;
#else
  return NULL;
#endif
}


void colvarproxy_python::py_check_obj(void *parent_obj,
                                      char const *obj_name,
                                      void *py_obj)
{
#if defined(COLVARS_PYTHON)
  if (py_obj == NULL) {
    std::string const obj_name_s(obj_name);
    PyObject *py_parent_obj =
      reinterpret_cast<PyObject *>((parent_obj != NULL) ? parent_obj :
                                   py_main);
    std::string parent_str = "";
    if (parent_obj != py_main) {
      PyObject *py_parent_name =
        PyObject_GetAttrString(py_parent_obj, "__name__");
      char const *py_parent_name_str =
        py_obj_to_str(reinterpret_cast<unsigned char *>(py_parent_name));
      parent_str = std::string(" within ")+std::string(py_parent_name_str);
      Py_XDECREF(py_parent_name);
    }
    cvm::error("Error: could not find Python attribute \""+
               std::string(obj_name_s)+"\""+parent_str+".\n", INPUT_ERROR);
  }
#endif
}


void *colvarproxy_python::py_new_array(size_t n)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_args = Py_BuildValue("(is)", n, "numpy.double");
  PyObject *py_v =
    PyObject_CallObject(reinterpret_cast<PyObject *>(py_np_zeros),
                        py_args);
  Py_DECREF(py_args);
  return reinterpret_cast<void *>(py_v);
#else
  return NULL;
#endif
}


void *colvarproxy_python::py_array_from_vector(std::vector<double> const &v,
                                               void *vector_obj)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_vector_obj = reinterpret_cast<PyObject *>(vector_obj);

  if (py_vector_obj == NULL) {
    py_vector_obj = reinterpret_cast<PyObject *>(py_new_array(v.size()));
  }
  // Copy into buffer
  Py_buffer py_vec_buf_view;
  PyObject_GetBuffer(py_vector_obj, &py_vec_buf_view, PyBUF_WRITABLE);
  double *const py_array = reinterpret_cast<double *>(py_vec_buf_view.buf);
  for (size_t i = 0; i < v.size(); i++) {
    py_array[i] = v[i];
  }
  PyBuffer_Release(&py_vec_buf_view);

  return reinterpret_cast<void *>(py_vector_obj);

#else

  return NULL;

#endif // #if defined(COLVARS_PYTHON)
}


void *colvarproxy_python::py_set_colvarvalue(colvarvalue const &value,
                                             void *value_obj)
{
#if defined(COLVARS_PYTHON)
  if (cvm::debug()) {
    cvm::log("colvarproxy_python::py_set_colvarvalue("+
             cvm::to_str(value)+").\n");
  }
  PyObject *py_value_obj = reinterpret_cast<PyObject *>(value_obj);

  if (py_value_obj == NULL) {
    if (value.type() == colvarvalue::type_scalar) {
      py_value_obj = PyFloat_FromDouble(0.0);
      if (cvm::debug()) {
        cvm::log("colvarproxy_python::py_set_colvarvalue(scalar).\n");
      }
    } else if (py_np_zeros != NULL) {
      // Create a new NumPy vector
      py_value_obj =
        reinterpret_cast<PyObject *>(py_new_array(value.size()));
      if (cvm::debug()) {
        cvm::log("colvarproxy_python::py_set_colvarvalue(new_vector).\n");
      }
    } else {
      // Create a new list
      py_value_obj = PyList_New(value.size());
      for (size_t ic = 0; ic < value.size(); ic++) {
        PyList_SetItem(py_value_obj, ic,
                       PyFloat_FromDouble(0.0)); // steals the ref
      }
      if (cvm::debug()) {
        cvm::log("colvarproxy_python::py_set_colvarvalue(list).\n");
      }
    }
  }

  if (value.type() == colvarvalue::type_scalar) {
    Py_DECREF(py_value_obj);
    py_value_obj = PyFloat_FromDouble(value.real_value);
  } else {
    // If py_value_obj is a NumPy array, the following will not be as fast as
    // copying into its underlying C array; however, retrieving the pointer
    // attribute from the PyObject could still be costlier for short vectors.
    for (size_t ic = 0; ic < value.size(); ic++) {
      PyObject *py_elem = PySequence_GetItem(py_value_obj, ic);
      if (py_elem) Py_DECREF(py_elem);
      py_elem = PyFloat_FromDouble(value[ic]);
      PySequence_SetItem(py_value_obj, ic, py_elem);
      Py_DECREF(py_elem);
    }
    if (cvm::debug()) {
      cvm::log("colvarproxy_python::py_set_colvarvalue(setting seq).\n");
    }
  }

  if (cvm::debug()) {
    cvm::log("colvarproxy_python::py_set_colvarvalue result = "+
             std::string(py_obj_to_str(
               reinterpret_cast<unsigned char *>(py_value_obj)))+".\n");
  }

  return reinterpret_cast<void *>(py_value_obj);

#else

  return NULL;

#endif
}


int colvarproxy_python::py_get_colvarvalue(void *value_obj,
                                           colvarvalue *value)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_value_obj = reinterpret_cast<PyObject *>(value_obj);
  size_t const n_dim = value->size();
  if (PyNumber_Check(py_value_obj) && (!PySequence_Check(py_value_obj))) {
    if (n_dim == 1) {
      *value = PyFloat_AsDouble(py_value_obj);
      return COLVARS_OK;
    }
  } else if (PySequence_Check(py_value_obj)) {
    if (PySequence_Size(py_value_obj) == int(n_dim)) {
      for (size_t ic = 0; ic < n_dim; ic++) {
        PyObject *py_num = PySequence_GetItem(py_value_obj, ic);
        (*value)[ic] = PyFloat_AsDouble(py_num);
        Py_DECREF(py_num);
      }
      return COLVARS_OK;
    }
  }

  PyErr_Print();
  return INPUT_ERROR;

#else

  return COLVARS_NOT_IMPLEMENTED;

#endif
}


void *colvarproxy_python::py_set_colvarvalues(
        std::vector<colvarvalue const *> const &values,
        void *values_obj)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_values_obj = reinterpret_cast<PyObject *>(values_obj);

  size_t const n_seq = values.size();

  if (py_values_obj == NULL) {
    py_values_obj = PyList_New(n_seq);
    for (size_t is = 0; is < n_seq; is++) {
      PyList_SetItem(py_values_obj, is, Py_None);
    }
  }

  for (size_t is = 0; is < n_seq; is++) {
    PyObject *py_elem =
      reinterpret_cast<PyObject *>(py_set_colvarvalue(*(values[is]), NULL));
    PyList_SetItem(py_values_obj, is, py_elem); // steals the ref
  }

  return reinterpret_cast<void *>(py_values_obj);

#else

  return NULL;

#endif
}


void *colvarproxy_python::py_set_colvar_gradients(
        std::vector<cvm::matrix2d<cvm::real> > const &grads,
        void *grads_obj)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_grads_obj = reinterpret_cast<PyObject *>(grads_obj);
  size_t const n_seq = grads.size();

  if (py_grads_obj == NULL) {
    py_grads_obj = PyList_New(n_seq);
    for (size_t is = 0; is < n_seq; is++) {
      PyList_SetItem(py_grads_obj, is, Py_None);
    }
  }

  for (size_t is = 0; is < n_seq; is++) {
    PyObject *py_elem =
      reinterpret_cast<PyObject *>(
        py_array_from_vector(grads[is].data_array(), NULL));
    PyList_SetItem(py_grads_obj, is, py_elem); // steals the ref
  }

  return reinterpret_cast<void *>(py_grads_obj);

#else

  return NULL;

#endif
}


int colvarproxy_python::py_get_colvar_gradients(
      void *grads_obj,
      std::vector<cvm::matrix2d<cvm::real> > *gradient)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_grads_obj = reinterpret_cast<PyObject *>(grads_obj);

  if (PyNumber_Check(py_grads_obj) && (!PySequence_Check(py_grads_obj))) {

    // a simple number is returned
    if ((gradient->size() == 1) && ((*gradient)[0].size() == 1)) {
      (*gradient)[0][0][0] = PyFloat_AsDouble(py_grads_obj);
      return COLVARS_OK;
    }

  } else {

    if (PySequence_Size(py_grads_obj) == int(gradient->size())) {

      for (size_t is = 0; is < gradient->size(); is++) {
        PyObject *py_elem = PySequence_GetItem(py_grads_obj, is);
        if (!PySequence_Check(py_elem)) {
          (*gradient)[is].data_array()[0] = PyFloat_AsDouble(py_elem);
          if (cvm::debug()) {
            cvm::log("Getting gradient of scalar component.\n");
            cvm::log("is = "+cvm::to_str(is)+".\n");
            cvm::log("grad = "+
                     cvm::to_str((*gradient)[is].data_array()[0], 22, 14)+
                     "\n");
          }
        } else {
          if (cvm::debug()) {
            cvm::log("Getting gradient of vector component.\n");
          }
          for (Py_ssize_t ic = 0; ic < PySequence_Size(py_elem); ic++) {
            // Copy matrix elements as a flattened vector
            // TODO optimize for vectors
            if (cvm::debug()) {
              cvm::log("is = "+cvm::to_str(is)+", ic = "+cvm::to_str(ic)+"\n");
            }
            PyObject *py_num = PySequence_GetItem(py_elem, ic);
            (*gradient)[is].data_array()[ic] = PyFloat_AsDouble(py_num);
            if (cvm::debug()) {
              cvm::log("grad = "+
                       cvm::to_str((*gradient)[is].data_array()[ic], 22, 14)+
                       "\n");
            }
            Py_DECREF(py_num);
          }
        }
        Py_DECREF(py_elem);
      }
      return COLVARS_OK;
    }
  }

  PyErr_Print();
  return INPUT_ERROR;

#else

  return COLVARS_NOT_IMPLEMENTED;

#endif
}


int colvarproxy_python::py_run_force_callback()
{
#if defined(COLVARS_PYTHON)
  if (py_main == NULL) init_py_pointers();
  PyGILState_STATE py_gil_state = PyGILState_Ensure();

  PyObject *py_func =
    reinterpret_cast<PyObject *>(py_get_function(py_main,
                                                 "calc_colvar_forces"));
  py_check_obj(py_main, "calc_colvar_forces", py_func);

  PyObject *py_step = PyLong_FromSize_t(cvm::step_absolute());
  PyObject *py_args = PyTuple_Pack(1, py_step);
  Py_DECREF(py_step);
  PyObject *py_result = PyObject_CallObject(py_func, py_args);
  Py_DECREF(py_args);
  if (py_result == NULL) {
    PyErr_Print();
    cvm::error("Error: in calling Python function calc_colvar_forces()\n",
               INPUT_ERROR);
    return cvm::get_error();
  } else {
    Py_DECREF(py_func);
  }
  // TODO add bias energy?

  PyGILState_Release(py_gil_state);
  return COLVARS_OK;

#else

  return COLVARS_NOT_IMPLEMENTED;

#endif
}


int colvarproxy_python::py_run_colvar_callback(void *func,
                                               void *cvc_values,
                                               void *colvar_value,
                                               void **out_obj)
{
#if defined(COLVARS_PYTHON)
  PyObject *py_func = reinterpret_cast<PyObject *>(func);
  PyObject *py_cvc_values = reinterpret_cast<PyObject *>(cvc_values);
  PyObject *py_colvar_value = reinterpret_cast<PyObject *>(colvar_value);
  PyObject *py_out_obj = reinterpret_cast<PyObject *>(*out_obj);

  PyObject *py_args =
    py_colvar_value ?
    PyTuple_Pack(3, py_cvc_values, py_colvar_value, py_out_obj) :
    PyTuple_Pack(2, py_cvc_values, py_out_obj);
  PyObject *py_result = PyObject_CallObject(py_func, py_args);
  Py_DECREF(py_args);

  if (py_result != NULL) {
    if (py_result != Py_None) {
      if (py_result != py_out_obj) {
        Py_XDECREF(py_out_obj);
        *out_obj = reinterpret_cast<void *>(py_result);
      }
    }
    return cvm::get_error();
  }

  PyErr_Print();
  return INPUT_ERROR;

#else

  return COLVARS_NOT_IMPLEMENTED;

#endif
}
