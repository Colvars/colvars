// -*- c++ -*-

#include <sstream>

#if defined(NAMD_TCL) || defined(VMDTCL)
#define COLVARS_TCL
#include <tcl.h>
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_tcl.h"
#include "colvaratoms.h"



colvarproxy_tcl::colvarproxy_tcl()
{
  tcl_interp_ = NULL;
}


colvarproxy_tcl::~colvarproxy_tcl()
{
}


void colvarproxy_tcl::init_tcl_pointers()
{
  cvm::error("Error: Tcl support is not available in this build.\n",
             COLVARS_NOT_IMPLEMENTED);
}


// TODO condense these via template when std::is_same is available

char const *colvarproxy_tcl::tcl_get_str(void *obj)
{
#if defined(COLVARS_TCL)
  return Tcl_GetString(reinterpret_cast<Tcl_Obj *>(obj));
#else
  return NULL;
#endif
}


int colvarproxy_tcl::tcl_get_int(void *obj, int *result)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *int_obj = reinterpret_cast<Tcl_Obj *>(obj);
  if (Tcl_GetIntFromObj(tcl_interp, int_obj, result) != TCL_OK) {
    return cvm::error("Error: could not get integer from Tcl object.\n",
                      INPUT_ERROR);
  }
  return COLVARS_OK;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_get_real(void *obj, cvm::real *result)
{
#if defined(COLVARS_TCL)
  double x = -1.0;
  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *int_obj = reinterpret_cast<Tcl_Obj *>(obj);
  if (Tcl_GetDoubleFromObj(tcl_interp, int_obj, &x) != TCL_OK) {
    cvm::error("Error: could not get real number from Tcl object.\n",
               INPUT_ERROR);
  }
  *result = x;
  return COLVARS_OK;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_get_int_vector(void *obj,
                                        std::vector<int> *vec)
{
#if defined(COLVARS_TCL)
  if (!vec) return BUG_ERROR;

  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *list_obj = reinterpret_cast<Tcl_Obj *>(obj);
  int list_length = 0;
  Tcl_Obj **obj_array = NULL;

  if (Tcl_ListObjGetElements(tcl_interp, list_obj, &list_length, &obj_array)
      == TCL_OK) {

    if (vec->size()) {
      if (vec->size() != list_length) {
        return cvm::error("Error: trying to copy a Tcl list into a vector "
                          "of different size.\n", INPUT_ERROR);
      }
    } else {
      vec->resize(list_length);
    }

    for (int i = 0; i < list_length; i++) {
      int x = 0;
      if (Tcl_GetIntFromObj(tcl_interp, obj_array[i], &x) !=
          TCL_OK) {
        return cvm::error("Error: could not read "+
                          cvm::to_str(i)+"-th element of Tcl list.\n",
                          INPUT_ERROR);
      } else {
        (*vec)[i] = x;
      }
    }
    return COLVARS_OK;
  } else {
    return cvm::error("Error: could not get values from Tcl list.\n",
                      INPUT_ERROR);
  }
  return COLVARS_OK;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_get_real_vector(void *obj,
                                         std::vector<cvm::real> *vec)
{
#if defined(COLVARS_TCL)
  if (!vec) return BUG_ERROR;

  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *list_obj = reinterpret_cast<Tcl_Obj *>(obj);
  int list_length = 0;
  Tcl_Obj **obj_array = NULL;

  if (Tcl_ListObjGetElements(tcl_interp, list_obj, &list_length, &obj_array)
      == TCL_OK) {

    if (vec->size()) {
      if (vec->size() != list_length) {
        return cvm::error("Error: trying to copy a Tcl list into a vector "
                          "of different size.\n", INPUT_ERROR);
      }
    } else {
      vec->resize(list_length);
    }

    for (int i = 0; i < list_length; i++) {
      double x = 0.0;
      if (Tcl_GetDoubleFromObj(tcl_interp, obj_array[i], &x) !=
          TCL_OK) {
        return cvm::error("Error: could not read "+
                          cvm::to_str(i)+"-th element of Tcl list.\n",
                          INPUT_ERROR);
      } else {
        (*vec)[i] = x;
      }
    }
    return COLVARS_OK;
  } else {
    return cvm::error("Error: could not get values from Tcl list.\n",
                      INPUT_ERROR);
  }
  return COLVARS_OK;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


void *colvarproxy_tcl::tcl_set_string(std::string const &s, void *obj)
{
#if defined(COLVARS_TCL)
  Tcl_Obj *tcl_obj = NULL;
  if (!obj) {
    tcl_obj = Tcl_NewStringObj(const_cast<char *>(s.c_str()), s.size());
  } else {
    tcl_obj = reinterpret_cast<Tcl_Obj *>(obj);
    Tcl_SetStringObj(tcl_obj, const_cast<char *>(s.c_str()), s.size());
  }
  return reinterpret_cast<void *>(tcl_obj);
#else
  return NULL;
#endif
}


void *colvarproxy_tcl::tcl_set_int(int const &x, void *obj)
{
#if defined(COLVARS_TCL)
  Tcl_Obj *tcl_obj = NULL;
  if (!obj) {
    tcl_obj = Tcl_NewIntObj(x);
  } else {
    tcl_obj = reinterpret_cast<Tcl_Obj *>(obj);
    Tcl_SetIntObj(tcl_obj, x);
  }
  return reinterpret_cast<void *>(tcl_obj);
#else
  return NULL;
#endif
}


void *colvarproxy_tcl::tcl_set_real(cvm::real const &x, void *obj)
{
#if defined(COLVARS_TCL)
  Tcl_Obj *tcl_obj = NULL;
  if (!obj) {
    tcl_obj = Tcl_NewDoubleObj(x);
  } else {
    tcl_obj = reinterpret_cast<Tcl_Obj *>(obj);
    Tcl_SetDoubleObj(tcl_obj, x);
  }
  return reinterpret_cast<void *>(tcl_obj);
#else
  return NULL;
#endif
}


void *colvarproxy_tcl::tcl_list_from_int_vec(std::vector<int> const &x)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *list = Tcl_NewListObj(0, NULL);
  for (size_t i = 0; i < x.size(); i++) {
    Tcl_ListObjAppendElement(interp, list, Tcl_NewIntObj(x[i]));
  }
  return reinterpret_cast<void *>(list);
#else
  return NULL;
#endif
}


void *colvarproxy_tcl::tcl_list_from_real_vec(std::vector<cvm::real> const &x)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *list = Tcl_NewListObj(0, NULL);
  for (size_t i = 0; i < x.size(); i++) {
    Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj((double) x[i]));
  }
  return reinterpret_cast<void *>(list);
#else
  return NULL;
#endif
}


void *colvarproxy_tcl::tcl_list_from_rvector_vec(std::vector<cvm::rvector> const &x)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *list = Tcl_NewListObj(0, NULL);
  for (size_t i = 0; i < x.size(); i++) {
    Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj((double) x[i].x));
    Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj((double) x[i].y));
    Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj((double) x[i].z));
  }
  return reinterpret_cast<void *>(list);
#else
  return NULL;
#endif
}


int colvarproxy_tcl::tcl_run_force_callback()
{
#if defined(COLVARS_TCL)
  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  std::string cmd = std::string("calc_colvar_forces ")
    + cvm::to_str(cvm::step_absolute());
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  if (err != TCL_OK) {
    cvm::log(std::string("Error while executing calc_colvar_forces:\n"));
    cvm::error(Tcl_GetStringResult(tcl_interp));
    return COLVARS_ERROR;
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_run_colvar_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         colvarvalue &value)
{
#if defined(COLVARS_TCL)

  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  size_t i;
  std::string cmd = std::string("calc_") + name;
  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") + (*(cvc_values[i])).to_simple_string() +
      std::string("}");
  }
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  const char *result = Tcl_GetStringResult(tcl_interp);
  if (err != TCL_OK) {
    return cvm::error(std::string("Error while executing ")
                      + cmd + std::string(":\n") +
                      std::string(Tcl_GetStringResult(tcl_interp)),
                      COLVARS_ERROR);
  }
  std::istringstream is(result);
  if (value.from_simple_string(is.str()) != COLVARS_OK) {
    cvm::log("Error parsing colvar value from script:");
    cvm::error(result);
    return COLVARS_ERROR;
  }
  return cvm::get_error();

#else

  return COLVARS_NOT_IMPLEMENTED;

#endif
}


int colvarproxy_tcl::tcl_run_colvar_gradient_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
#if defined(COLVARS_TCL)

  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  size_t i;
  std::string cmd = std::string("calc_") + name + "_gradient";
  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") + (*(cvc_values[i])).to_simple_string() +
      std::string("}");
  }
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  if (err != TCL_OK) {
    return cvm::error(std::string("Error while executing ")
                      + cmd + std::string(":\n") +
                      std::string(Tcl_GetStringResult(tcl_interp)),
                      COLVARS_ERROR);
  }
  Tcl_Obj **list;
  int n;
  Tcl_ListObjGetElements(tcl_interp, Tcl_GetObjResult(tcl_interp),
                         &n, &list);
  if (n != int(gradient.size())) {
    cvm::error("Error parsing list of gradient values from script: found "
               + cvm::to_str(n) + " values instead of " +
               cvm::to_str(gradient.size()));
    return COLVARS_ERROR;
  }
  for (i = 0; i < gradient.size(); i++) {
    std::istringstream is(Tcl_GetString(list[i]));
    if (gradient[i].from_simple_string(is.str()) != COLVARS_OK) {
      cvm::log("Gradient matrix size: " + cvm::to_str(gradient[i].size()));
      cvm::log("Gradient string: " + cvm::to_str(Tcl_GetString(list[i])));
      cvm::error("Error parsing gradient value from script", COLVARS_ERROR);
      return COLVARS_ERROR;
    }
  }

  return cvm::get_error();

#else

  return COLVARS_NOT_IMPLEMENTED;

#endif
}
