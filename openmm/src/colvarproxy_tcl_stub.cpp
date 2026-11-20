#include <string>
#include <vector>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarvalue.h"
#include "colvarproxy_tcl.h"

colvarproxy_tcl::colvarproxy_tcl()
    : tcl_interp_(nullptr)
{
}

colvarproxy_tcl::~colvarproxy_tcl() = default;

char const* colvarproxy_tcl::tcl_get_str(void*)
{
    return "";
}

int colvarproxy_tcl::tcl_run_script(std::string const&)
{
    return cvm::error_static("Error: Tcl scripting is not supported by the OpenMM Colvars interface.\n",
                             COLVARS_NOT_IMPLEMENTED);
}

int colvarproxy_tcl::tcl_run_file(std::string const&)
{
    return cvm::error_static("Error: Tcl scripting is not supported by the OpenMM Colvars interface.\n",
                             COLVARS_NOT_IMPLEMENTED);
}

int colvarproxy_tcl::tcl_run_force_callback()
{
    return cvm::error_static("Error: Tcl force callbacks are not supported by the OpenMM Colvars interface.\n",
                             COLVARS_NOT_IMPLEMENTED);
}

int colvarproxy_tcl::tcl_run_colvar_callback(std::string const&,
                                             std::vector<const colvarvalue*> const&,
                                             colvarvalue&)
{
    return cvm::error_static("Error: Tcl colvar callbacks are not supported by the OpenMM Colvars interface.\n",
                             COLVARS_NOT_IMPLEMENTED);
}

int colvarproxy_tcl::tcl_run_colvar_gradient_callback(
    std::string const&,
    std::vector<const colvarvalue*> const&,
    std::vector<cvm::matrix2d<cvm::real>>&)
{
    return cvm::error_static("Error: Tcl colvar gradient callbacks are not supported by the OpenMM Colvars interface.\n",
                             COLVARS_NOT_IMPLEMENTED);
}

void colvarproxy_tcl::init_tcl_pointers()
{
}
