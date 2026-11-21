#include <iostream>
#include <string>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvarscript.h"
#include "colvar.h"
#include "colvarcomp.h"

namespace {
class tabf_proxy : public colvarproxy_stub {
public:
  bool fail_boundaries = false;
  int boundary_queries = 0;
  int sends = 0;
  void start_simulation() { b_simulation_running = true; }
  bool reject_send = false;
  bool reject_frequency = false;

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
  int request_tabf_energy_freq(int const) override {
    if (reject_frequency) {
      cvmodule->error("TABF requires a different energy schedule.");
      return COLVARS_INPUT_ERROR;
    }
    return colvarproxy_system::request_tabf_energy_freq(1);
  }
  int send_tabf_alpha() override {
    ++sends;
    return reject_send ? COLVARS_NOT_IMPLEMENTED : COLVARS_OK;
  }
};

bool reported_once(colvarproxy &proxy, std::string const &message) {
  auto const &errors = proxy.get_error_msgs();
  auto const first = errors.find(message);
  return first != std::string::npos && errors.find(message, first + 1) == std::string::npos;
}

bool check_estimators(bool mixed, std::string const &options, std::string const &error = "") {
  tabf_proxy proxy;
  proxy.set_unit_system("real", false);
  std::string conf = "colvar {\n name a\n width 0.01\n tabfAlpha {}\n}\n";
  if (mixed) {
    conf += "colvar {\n name d\n distance {\n group1 { atomNumbers 1 }\n group2 { atomNumbers 2 }\n }\n"
            " upperBoundary 10\n extendedLagrangian on\n extendedTemp 300\n extendedFluctuation 1\n}\n";
  }
  conf += std::string("abf {\n colvars a") + (mixed ? " d" : "") +
    "\n applyBias off\n updateBias off\n" + options + "\n}\n";
  int const rc = proxy.cvmodule->read_config_string(conf);
  return error.empty() ? rc == COLVARS_OK : ((rc & COLVARS_INPUT_ERROR) && reported_once(proxy, error));
}

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
      ok = ok && cv->lower_boundary.real_value == lower && cv->upper_boundary.real_value == upper;
      if (components == "tabfAlpha {}") ok = ok && proxy.boundary_queries == 1;
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
  failures += !check_case("within tolerance uses exact limits", "lowerBoundary 0.9700000000005\n");
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
  failures += !check_estimators(false, "");
  failures += !check_estimators(false, "UIestimator on", "UIestimator cannot be used");
  failures += !check_estimators(false, "CZARestimator on", "CZARestimator cannot be used");
  failures += !check_estimators(true, "", "CZARestimator cannot be used"); // eABF defaults CZAR on
  failures += !check_estimators(true, "CZARestimator off\nUIestimator off");

  // Base proxy methods only return status; the owning caller reports it.
  {
    colvarproxy_stub proxy;
    cvm::real alpha, lower, upper;
    bool ok = proxy.get_tabf_alpha(&alpha) == COLVARS_NOT_IMPLEMENTED &&
      proxy.get_tabf_alpha_boundaries(&lower, &upper) == COLVARS_NOT_IMPLEMENTED &&
      proxy.get_dE_dtabf_alpha(&alpha) == COLVARS_NOT_IMPLEMENTED &&
      proxy.send_tabf_alpha() == COLVARS_NOT_IMPLEMENTED &&
      proxy.request_tabf_energy_freq(1) == COLVARS_OK && proxy.get_error_msgs().empty();
    int const rc = proxy.cvmodule->read_config_string("colvar {\n name a\n tabfAlpha {}\n}\n");
    ok = ok && (rc & COLVARS_NOT_IMPLEMENTED) && reported_once(proxy, "requires a NAMD build with TABF support");
    std::cout << (ok ? "PASS " : "FAIL ") << "unsupported proxy\n";
    failures += !ok;
  }
  {
    tabf_proxy proxy;
    int rc = proxy.cvmodule->read_config_string("colvar {\n name a\n tabfAlpha {}\n"
                                                " extendedLagrangian on\n extendedMass 1\n}\n");
    colvar *cv = proxy.cvmodule->colvar_by_name("a");
    bool ok = rc == COLVARS_OK && cv;
    if (ok) {
      rc = proxy.cvmodule->calc(); // proxy supplies alpha, but not dE/dalpha
      ok = (rc & COLVARS_NOT_IMPLEMENTED) && reported_once(proxy, "get_dE_dtabf_alpha");
      proxy.cvmodule->clear_error();
      proxy.start_simulation();
      unsigned char cmd[] = "cv", sub[] = "colvar", name[] = "a", update[] = "update";
      unsigned char *args[] = {cmd, sub, name, update};
      cvm::real const before = cv->value().real_value;
      rc = proxy.script->run(4, args);
      ok = ok && rc != COLVARS_OK && reported_once(proxy, "get_dE_dtabf_alpha") &&
        cv->value().real_value == before && proxy.end_of_step() == COLVARS_OK && proxy.sends == 0;
      proxy.cvmodule->clear_error();
      proxy.reject_send = true;
      rc = cv->set_state_params("x 1\nextended_x 1\nextended_v 0\n");
      ok = ok && (rc & COLVARS_NOT_IMPLEMENTED) && reported_once(proxy, "send_tabf_alpha");
      proxy.cvmodule->clear_error();
      proxy.reject_send = false;
      ok = ok && proxy.end_of_step() == COLVARS_OK && proxy.sends == 2; // retry
      ok = ok && proxy.end_of_step() == COLVARS_OK && proxy.sends == 2; // clean cache
      proxy.set_tabf_alpha(1.2);
      proxy.reject_send = true;
      ok = ok && (proxy.end_of_step() & COLVARS_NOT_IMPLEMENTED) && reported_once(proxy, "send_tabf_alpha");
    }
    std::cout << (ok ? "PASS " : "FAIL ") << "runtime read, restart and deferred send\n";
    failures += !ok;
  }
  // Keep a specific back-end diagnostic and its returned status without a second message.
  {
    tabf_proxy proxy;
    proxy.reject_frequency = true;
    int const rc = proxy.cvmodule->read_config_string("colvar {\n name a\n tabfAlpha {}\n}\n");
    bool const ok = (rc & COLVARS_INPUT_ERROR) && reported_once(proxy, "different energy schedule") &&
      proxy.get_error_msgs().find("request_tabf_energy_freq") == std::string::npos;
    std::cout << (ok ? "PASS " : "FAIL ") << "back-end frequency error\n";
    failures += !ok;
  }
  // Test the failing query directly: no intrinsic features may be registered,
  // even when the backend returns an error without writing either output.
  {
    tabf_proxy proxy;
    proxy.fail_boundaries = true;
    colvar::tabf_alpha component;
    proxy.cvmodule->error("Unrelated error already recorded.");
    int const result = component.init("");
    bool const ok = result == COLVARS_NOT_IMPLEMENTED && proxy.boundary_queries == 1 &&
      reported_once(proxy, "get_tabf_alpha_boundaries") &&
      !component.is_enabled(colvardeps::f_cvc_lower_boundary) &&
      !component.is_enabled(colvardeps::f_cvc_upper_boundary);
    std::cout << (ok ? "PASS " : "FAIL ") << "backend query failure\n";
    failures += !ok;
  }
  return failures ? 1 : 0;
}
