#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvar.h"
#include "colvarbias.h"

namespace {
void require(bool condition, std::string const &message)
{
  if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, std::string const &message)
{
  require(std::fabs(actual - expected) < 1.0e-10 * (1.0 + std::fabs(expected)),
          message + ": actual=" + std::to_string(actual) + ", expected=" + std::to_string(expected));
}

class meta_proxy : public colvarproxy_stub {
public:
  double energy = 0.0;
  void add_energy(cvm::real value) override { energy += value; }
  void set_distance(double distance)
  {
    auto &positions = *modify_atom_positions();
    positions[0] = cvm::rvector(0.0, 0.0, 0.0);
    positions[1] = cvm::rvector(distance, 0.0, 0.0);
  }
  void reset_buffers()
  {
    energy = 0.0;
    for (auto &force : *modify_atom_applied_forces()) force.reset();
  }
  double force0() const { return (*get_atom_applied_forces())[0].x; }
};

std::string configuration(bool grids)
{
  return std::string("colvar {\n name d\n width 0.05\n lowerBoundary 0.0\n upperBoundary 3.0\n"
                     " distance {\n group1 { atomNumbers 1 }\n group2 { atomNumbers 2 }\n }\n}\n"
                     "metadynamics {\n name mtd\n colvars d\n hillWeight 0.5\n hillWidth 2.0\n"
                     " newHillFrequency 2\n useGrids ") + (grids ? "on" : "off") + "\n}\n";
}

void test_trial_matches_commit(bool grids)
{
  meta_proxy proxy;
  require(proxy.cvmodule->read_config_string(configuration(grids)) == COLVARS_OK, "configuration failed");
  colvarbias *bias = proxy.cvmodule->bias_by_name("mtd");
  require(bias != nullptr, "metadynamics bias not found");
  require(bias->trial_mode() == colvarbias::trial_evaluate, "metadynamics is not trial-evaluated");

  // Committed steps 0-4 deposit hills at steps 0, 2 and 4 (and project them on the grid)
  double const history[] = {1.00, 1.03, 1.06, 1.09, 1.12};
  for (double distance : history) {
    proxy.set_distance(distance);
    proxy.reset_buffers();
    require(proxy.cvmodule->calc() == COLVARS_OK, "committed step failed");
    proxy.cvmodule->it++;
  }

  std::string state_before;
  require(bias->write_state_string(state_before) == COLVARS_OK, "state output failed");
  double const committed_energy = bias->get_energy();
  cvm::step_number const step = proxy.cvmodule->it;

  // Step 5 deposits no hill and updates no grid: trial and commit must agree there
  double const distance = 1.08;
  double trial_energy = 0.0, trial_force = 0.0;
  for (double trial_distance : {distance, 1.20, distance}) {
    proxy.set_distance(trial_distance);
    proxy.reset_buffers();
    require(proxy.cvmodule->calc_trial() == COLVARS_OK, "trial failed");
    require(!proxy.cvmodule->in_trial(), "trial flag leaked");
    std::string state_after;
    require(bias->write_state_string(state_after) == COLVARS_OK, "state output failed");
    require(state_after == state_before, "trial changed the hills or metadynamics state");
    near(bias->get_energy(), committed_energy, "trial did not restore the committed bias energy");
    require(proxy.cvmodule->it == step, "trial changed the step");
    trial_energy = proxy.energy;
    trial_force = proxy.force0();
  }
  require(trial_energy > 1.0e-3, "no hill contributes: the test is trivial");
  require(std::fabs(trial_force) > 1.0e-6, "no hill force: the test is trivial");

  proxy.set_distance(distance);
  proxy.reset_buffers();
  require(proxy.cvmodule->calc() == COLVARS_OK, "committed comparison step failed");
  near(trial_energy, proxy.energy, "trial energy differs from the committed energy");
  near(trial_force, proxy.force0(), "trial force differs from the committed force");
}
}

int main()
{
  bool failed = false;
  for (bool grids : {false, true}) {
    char const *name = grids ? "useGrids on" : "useGrids off";
    try {
      test_trial_matches_commit(grids);
      std::cout << "PASS: metadynamics trial (" << name << ")\n";
    } catch (std::exception const &error) {
      failed = true;
      std::cerr << "FAIL: metadynamics trial (" << name << "): " << error.what() << '\n';
    }
  }
  return failed ? 1 : 0;
}
