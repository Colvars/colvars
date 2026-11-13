#include <cmath>
#include <iostream>
#include <string>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvar.h"
#include "colvarcomp.h"

namespace {
class tabf_proxy : public colvarproxy_stub {
public:
  bool fail_boundaries = false;
  int boundary_queries = 0;

  int get_tabf_alpha(cvm::real *alpha) override {
    *alpha = 1.0;
    return COLVARS_OK;
  }
  int get_tabf_alpha_boundaries(cvm::real *lower, cvm::real *upper) override {
    ++boundary_queries;
    if (fail_boundaries) return COLVARS_NOT_IMPLEMENTED;
    *lower = 0.97;
    *upper = 2.03;
    return COLVARS_OK;
  }
  int request_tabf_energy_freq(int const) override { return COLVARS_OK; }
};

bool close(cvm::real a, cvm::real b) { return std::fabs(a - b) < 1.0e-12; }

bool check_case(std::string const &name, std::string const &options,
                std::string const &expected_error = "", bool hard = true,
                std::string const &components = "tabfAlpha {}",
                bool intrinsic = true, cvm::real lower = 0.97,
                cvm::real upper = 2.03) {
  tabf_proxy proxy;
  proxy.set_unit_system("real", false);
  int const result = proxy.cvmodule->read_config_string(
    "colvar {\n name alpha\n width 0.01\n" + options + "\n" + components + "\n}\n");
  bool const success = expected_error.empty();
  bool ok = (result == COLVARS_OK) == success;
  if (!success) {
    ok = ok && (result & COLVARS_INPUT_ERROR) &&
      proxy.get_error_msgs().find(expected_error) != std::string::npos;
  }
  if (success && result == COLVARS_OK) {
    colvar *cv = proxy.cvmodule->colvar_by_name("alpha");
    ok = ok && cv &&
      cv->is_enabled(colvardeps::f_cv_lower_boundary) == intrinsic &&
      cv->is_enabled(colvardeps::f_cv_upper_boundary) == intrinsic &&
      cv->is_enabled(colvardeps::f_cv_hard_lower_boundary) == hard &&
      cv->is_enabled(colvardeps::f_cv_hard_upper_boundary) == hard;
    if (cv && intrinsic) {
      ok = ok && close(cv->lower_boundary.real_value, lower) &&
        close(cv->upper_boundary.real_value, upper);
    }
  }
  std::cout << (ok ? "PASS " : "FAIL ") << name << '\n';
  return ok;
}
}

int main() {
  int failures = 0;
  std::string const on = "hardLowerBoundary on\nhardUpperBoundary on\n";
  failures += !check_case("both omitted, hard on", on);
  failures += !check_case("lower omitted", on + "upperBoundary 2.03\n");
  failures += !check_case("upper omitted", on + "lowerBoundary 0.97\n");
  failures += !check_case("intrinsic hard defaults", "");
  failures += !check_case("explicit matching", "lowerBoundary 0.97\nupperBoundary 2.03\n");
  failures += !check_case("conflicting lower", "lowerBoundary 0.98\n",
                          "tabfAlpha lowerBoundary 0.98 conflicts with back-end TABF alpha minimum 0.97");
  failures += !check_case("conflicting upper", "upperBoundary 2.02\n",
                          "tabfAlpha upperBoundary 2.02 conflicts with back-end TABF alpha maximum 2.03");
  failures += !check_case("explicit soft boundaries", "hardLowerBoundary off\nhardUpperBoundary off\n", "", false);
  failures += !check_case("scaled component has no alpha defaults", "", "", false,
                          "tabfAlpha {\n componentCoeff 2\n}", false);
  failures += !check_case("nonlinear component has no alpha defaults", "", "", false,
                          "tabfAlpha {\n componentExp 2\n}", false);
  failures += !check_case("multiple components have no alpha defaults", "", "", false,
                          "tabfAlpha {}\ntabfAlpha {}", false);
  failures += !check_case("non-TABF bounded component", "", "", true,
                          "angle {\n group1 { atomNumbers 1 }\n group2 { atomNumbers 2 }\n group3 { atomNumbers 3 }\n}",
                          true, 0.0, 180.0);
  // Test the failing query directly: no intrinsic features may be registered,
  // even when the backend returns an error without writing either output.
  {
    tabf_proxy proxy;
    proxy.fail_boundaries = true;
    colvar::tabf_alpha component;
    int const result = component.init("");
    bool const ok = result == COLVARS_NOT_IMPLEMENTED && proxy.boundary_queries == 1 &&
      !component.is_enabled(colvardeps::f_cvc_lower_boundary) &&
      !component.is_enabled(colvardeps::f_cvc_upper_boundary);
    std::cout << (ok ? "PASS " : "FAIL ") << "backend query failure\n";
    failures += !ok;
  }
  return failures ? 1 : 0;
}
