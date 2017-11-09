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
  _tcl_interp = NULL;
}


colvarproxy_tcl::~colvarproxy_tcl()
{
}


void colvarproxy_tcl::init_tcl_pointers()
{
  cvm::error("Error: Tcl support is currently unavailable "
             "outside NAMD or VMD.\n", COLVARS_NOT_IMPLEMENTED);
}


char const *colvarproxy_tcl::tcl_get_str(void *obj)
{
#if defined(COLVARS_TCL)
  return Tcl_GetString(reinterpret_cast<Tcl_Obj *>(obj));
#else
  return NULL;
#endif
}


int colvarproxy_tcl::tcl_get_int(void *obj)
{
  int result = -1;
#if defined(COLVARS_TCL)
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *int_obj = reinterpret_cast<Tcl_Obj *>(obj);
  if (Tcl_GetIntFromObj(interp, int_obj, &result) != TCL_OK) {
    cvm::error("Error: could not get integer from Tcl object.\n",
               INPUT_ERROR);
  }
#endif
  return result;
}


std::vector<int> colvarproxy_tcl::tcl_get_int_vector(void *obj)
{
  std::vector<int> result;
#if defined(COLVARS_TCL)
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *list_obj = reinterpret_cast<Tcl_Obj *>(obj);
  int list_length = 0;
  Tcl_Obj **int_obj_array = NULL;
  if (Tcl_ListObjGetElements(interp, list_obj, &list_length, &int_obj_array)
      == TCL_OK) {
    result.resize(list_length, 0);
    for (int i = 0; i < list_length; i++) {
      if (Tcl_GetIntFromObj(interp, int_obj_array[i], &(result[i])) != TCL_OK) {
        cvm::error("Error: could not get integer from index "+
                   cvm::to_str(i)+" of Tcl list.\n",
                   INPUT_ERROR);
        break;
      }
    }
  } else {
    cvm::error("Error: could not get values from Tcl list.\n",
               INPUT_ERROR);
  }
#endif
  return result;
}


void *colvarproxy_tcl::tcl_set_string(std::string const &s)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *tcl_x = Tcl_NewStringObj(const_cast<char *>(s.c_str()), s.size());
  return tcl_x;
#else
  return NULL;
#endif
}


void *colvarproxy_tcl::tcl_set_real(cvm::real x)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *tcl_x = Tcl_NewDoubleObj((double) x);
  return tcl_x;
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
  return list;
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
  return list;
#else
  return NULL;
#endif
}


int colvarproxy_tcl::tcl_get_real(void *obj, cvm::real *x)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *interp = reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  Tcl_Obj *tcl_obj = reinterpret_cast<Tcl_Obj *>(obj);
  double y = 0.0;
  if (obj) {
    Tcl_GetDoubleFromObj(interp, tcl_obj, &y);
    *x = y;
  } else {
    return INPUT_ERROR;
  }
  return COLVARS_OK;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_run_force_callback()
{
#if defined(COLVARS_TCL)
  Tcl_Interp *const tcl_interp = reinterpret_cast<Tcl_Interp *>(_tcl_interp);
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

  Tcl_Interp *const tcl_interp = reinterpret_cast<Tcl_Interp *>(_tcl_interp);
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

  Tcl_Interp *const tcl_interp = reinterpret_cast<Tcl_Interp *>(_tcl_interp);
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


int colvarproxy_tcl::tcl_add_atoms_callback(void *atom_list,
                                            cvm::atom_group *atoms)
{
#if defined(COLVARS_TCL)
  return atoms->add_atom_numbers(tcl_get_int_vector(atom_list));
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}
