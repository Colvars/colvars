%module openmm_colvars
%pythonbegin %{
import openmm
%}
%{
#include "ColvarsForce.h"
#include "openmm/Context.h"
#include <exception>
%}
%include "std_string.i"
%include "std_vector.i"
%import(module="openmm.openmm") "openmm_colvars_openmm_import.i"
namespace std { %template(vector_double) vector<double>; }
%exception {
    try { $action }
    catch (const std::exception& exception) {
        PyErr_SetString(PyExc_RuntimeError, exception.what());
        SWIG_fail;
    }
    catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown exception in OpenMM Colvars");
        SWIG_fail;
    }
}
%include "ColvarsForce.h"
