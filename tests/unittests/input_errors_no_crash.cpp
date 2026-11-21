// Check that invalid inputs are reported as errors without crashing when the
// proxy's error() function returns, as it does in VMD and in the stub proxy
// (MD engines usually stop the run inside error(), which hides such bugs)

#include <iostream>
#include <string>
#include <vector>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"

namespace {

std::string const colvars_def = R"(
colvar {
  name a
  distanceDir {
    group1 {
      atomNumbers 1 2
    }
    group2 {
      atomNumbers 3 4
    }
  }
}

colvar {
  name b
  distance {
    group1 {
      atomNumbers 5 6
    }
    group2 {
      atomNumbers 7 8
    }
  }
}
)";

struct test_case {
  std::string description;
  std::string config;
};

std::vector<test_case> const test_cases = {
  { "targetCenters with the wrong number of values",
    R"(
harmonic {
  colvars a
  centers (1.0, 0.0, 0.0)
  forceConstant 1.0
  targetCenters 1.0 2.0
  targetNumSteps 2
}
)" },
  { "maxForce with the wrong number of values",
    R"(
abf {
  colvars b a
  maxForce 1.0
}
)" },
};

} // namespace


int main(int argc, char *argv[]) {

  int retval = 0;

  for (auto const &tc : test_cases) {

    colvarproxy_stub *proxy = new colvarproxy_stub();
    proxy->set_unit_system("real", false);
    for (int ai = 0; ai < 104; ai++) {
      proxy->init_atom(ai+1);
    }

    int err = proxy->cvmodule->read_config_string(colvars_def);
    if (err != COLVARS_OK) {
      std::cerr << "Error: could not define the variables for test \"" << tc.description
                << "\".\n";
      retval = 1;
    } else {
      err = proxy->cvmodule->read_config_string(tc.config);
      if (err == COLVARS_OK) {
        std::cerr << "Error: invalid input was not detected: " << tc.description << ".\n";
        retval = 1;
      }
    }

    delete proxy;
  }

  return retval;
}
