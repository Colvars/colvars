#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvarcomp_coordnums.h"

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

class pairlist_proxy : public colvarproxy_stub {
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
};

struct pairlist_probe : colvar::coordnum {
  static bool entry(colvar::coordnum const &cvc) { return (cvc.*&pairlist_probe::pairlist)[0]; }
};

std::string configuration(bool self)
{
  std::string const component = self ?
    " selfCoordNum {\n group1 { atomNumbers 1 2 }\n" :
    " coordNum {\n group1 { atomNumbers 1 }\n group2 { atomNumbers 2 }\n";
  return "colvar {\n name cn\n width 1\n" + component +
    " cutoff 1\n tolerance 0.1\n pairListFrequency 10\n }\n}\n"
    "harmonic {\n colvars cn\n centers 0\n forceConstant 2\n}\n";
}

void test_contracted_trial(bool self)
{
  double trial_value = 0.0, trial_energy = 0.0, trial_force = 0.0;
  double const distance = 0.8;
  {
    pairlist_proxy proxy;
    require(proxy.cvmodule->read_config_string(configuration(self)) == COLVARS_OK, "configuration failed");
    proxy.set_distance(4.0);
    require(proxy.cvmodule->calc() == COLVARS_OK, "committed rebuild failed");
    auto *cv = proxy.cvmodule->colvar_by_name("cn");
    auto const &cvc = dynamic_cast<colvar::coordnum const &>(*cv->get_cvcs()[0]);
    require(!pairlist_probe::entry(cvc), "far pair was not excluded by committed rebuild");
    proxy.cvmodule->it = 1;
    std::string const state = cv->get_state_params();
    for (double trial_distance : {distance, 4.0, distance}) {
      proxy.set_distance(trial_distance);
      proxy.reset_buffers();
      require(proxy.cvmodule->calc_trial() == COLVARS_OK, "contracted trial failed");
      require(!pairlist_probe::entry(cvc), "trial mutated persistent pairlist");
      require(cv->get_state_params() == state, "trial mutated committed CV history");
      require(proxy.cvmodule->it == 1 && !proxy.cvmodule->in_trial(), "trial changed step or leaked trial flag");
    }
    trial_value = cvc.value().real_value;
    trial_energy = proxy.energy;
    trial_force = (*proxy.get_atom_applied_forces())[0].x;
  }
  // A direct committed evaluation rebuilds its list at the contracted coordinates.
  {
    pairlist_proxy proxy;
    require(proxy.cvmodule->read_config_string(configuration(self)) == COLVARS_OK, "direct configuration failed");
    proxy.set_distance(distance);
    require(proxy.cvmodule->calc() == COLVARS_OK, "direct evaluation failed");
    auto *cv = proxy.cvmodule->colvar_by_name("cn");
    double const s = 1.0 / (1.0 + std::pow(distance, 6));
    double const value = (s - 0.1) / 0.9;
    double const derivative = -6.0 * std::pow(distance, 5) * s * s / 0.9;
    near(cv->value().real_value, value, "direct tolerance-corrected CV differs from analytic value");
    near(proxy.energy, value * value, "direct energy differs from analytic energy");
    near((*proxy.get_atom_applied_forces())[0].x, 2.0 * value * derivative, "direct force differs from analytic force");
    near(trial_value, value, "contracted trial omitted an entering pair");
    near(trial_energy, proxy.energy, "trial energy differs from direct rebuild");
    near(trial_force, (*proxy.get_atom_applied_forces())[0].x, "trial force differs from direct rebuild");
  }
}
}

int main()
{
  bool failed = false;
  for (bool self : {false, true}) {
    char const *name = self ? "selfCoordNum" : "coordNum";
    try {
      test_contracted_trial(self);
      std::cout << "PASS: " << name << " contracted trial and persistent pairlist/history\n";
    } catch (std::exception const &error) {
      failed = true;
      std::cerr << "FAIL: " << name << ": " << error.what() << '\n';
    }
  }
  return failed ? 1 : 0;
}
