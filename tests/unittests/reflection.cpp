#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvar.h"

// Runtime integer inputs keep the test independent of floating classification
// and constant folding, even when this entire translation unit uses fast math.
double from_bits(std::uint64_t bits) {
  volatile std::uint64_t runtime = bits;
  bits = runtime;
  double value;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

int main() {
  colvarproxy_stub proxy;
  colvar cv(proxy.cvmodule);
  cv.set_enabled(colvardeps::f_cv_scalar, true);
  cv.set_enabled(colvardeps::f_cv_lower_boundary, true);
  cv.set_enabled(colvardeps::f_cv_upper_boundary, true);
  cv.set_enabled(colvardeps::f_cv_extended_Lagrangian, true);
  int failures = 0;
  auto check = [&](char const *name, double value, double lower, double upper,
                   bool lo, bool hi, double expected, bool bounced, bool odd,
                   char const *error = nullptr) {
    cv.lower_boundary = colvarvalue(lower);
    cv.upper_boundary = colvarvalue(upper);
    // Directly exercise the production method, including invalid intervals
    // that the configuration parser would normally reject earlier.
    cv.set_enabled(colvardeps::f_cv_reflecting_lower_boundary, lo);
    cv.set_enabled(colvardeps::f_cv_reflecting_upper_boundary, hi);
    colvarvalue position(value);
    bool reflected = false, parity = false;
    proxy.cvmodule->clear_error();
    proxy.clear_error_msgs();
    int result = cv.reflect_extended_coordinate(position, reflected, parity);
    bool ok = error ? ((result & COLVARS_INPUT_ERROR) &&
      proxy.get_error_msgs().find(error) != std::string::npos) :
      (result == COLVARS_OK && std::fabs(position.real_value - expected) < 1.e-12 &&
       reflected == bounced && parity == odd);
    std::cout << (ok ? "PASS " : "FAIL ") << name << '\n';
    failures += !ok;
  };
  check("inside", 0.5, 0, 1, true, true, 0.5, false, false);
  check("lower endpoint", 0, 0, 1, true, true, 0, false, false);
  check("upper endpoint", 1, 0, 1, true, true, 1, false, false);
  check("lower only", -0.25, 0, 1, true, false, 0.25, true, true);
  check("upper only", 1.25, 0, 1, false, true, 0.75, true, true);
  check("no reflection", -2, 0, 1, false, false, -2, false, false);
  check("shifted interval", 13, 10, 12, true, true, 11, true, true);
  check("one crossing", 1.25, 0, 1, true, true, 0.75, true, true);
  check("two crossings", 2.25, 0, 1, true, true, 0.25, true, false);
  check("many positive", 101.25, 0, 1, true, true, 0.75, true, true);
  check("many negative", -100.25, 0, 1, true, true, 0.25, true, true);
  check("positive lower landing", 4, 0, 1, true, true, 0, true, true);
  check("positive upper landing", 5, 0, 1, true, true, 1, true, false);
  check("negative lower landing", -4, 0, 1, true, true, 0, true, false);
  check("negative upper landing", -3, 0, 1, true, true, 1, true, true);
  for (std::uint64_t bits : {UINT64_C(0x7ff8000000001234),
                            UINT64_C(0x7ff0000000000000),
                            UINT64_C(0xfff0000000000000)}) {
    double invalid = from_bits(bits);
    check("nonfinite input", invalid, 0, 1, true, true, 0, false, false, "is non-finite");
    check("nonfinite lower", 0.5, invalid, 1, true, true, 0, false, false, "finite, positive interval");
    check("nonfinite upper", 0.5, 0, invalid, true, true, 0, false, false, "finite, positive interval");
  }
  check("zero interval", 2, 1, 1, true, true, 0, false, false, "finite, positive interval");
  check("reversed interval", 2, 1, 0, true, true, 0, false, false, "finite, positive interval");
  double const max = from_bits(UINT64_C(0x7fefffffffffffff));
  check("width overflow", 0, -max, max, true, true, 0, false, false, "finite, positive interval");
  check("period overflow", 0, 0, max, true, true, 0, false, false, "finite, positive interval");
  check("shift overflow", max, -max, -max * 0.75, true, true, 0, false, false, "overflow while applying reflecting boundaries");
  check("single overflow", -max, max * 0.75, max, true, false, 0, false, false, "overflow while applying reflecting boundary");
  return failures ? 1 : 0;
}
