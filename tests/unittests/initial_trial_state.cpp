#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvar.h"

namespace {
void require(bool condition, std::string const &message)
{
  if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, std::string const &message)
{
  require(std::fabs(actual - expected) < 1.0e-10 * (1.0 + std::fabs(expected)), message);
}

class trial_proxy : public colvarproxy_stub {
public:
  double energy = 0.0;
  double lambda = 0.37;
  int random_calls = 0;
  int derivative_calls = 0;
  int writeback_calls = 0;
  bool throw_energy = false;
  bool throw_lambda = false;

  trial_proxy()
  {
    b_simulation_running = true;
    set_unit_system("real", false);
    set_target_temperature(300.0);
    set_integration_timestep(1.0);
    cvm::rotation::monitor_crossings = false;
  }
  void add_energy(cvm::real value) override
  {
    if (throw_energy) throw std::runtime_error("injected trial failure");
    energy += value;
  }
  cvm::real rand_gaussian() override { ++random_calls; return 0.25; }
  int get_alch_lambda(cvm::real *value) override
  {
    if (throw_lambda) throw std::runtime_error("injected preparation failure");
    *value = lambda;
    return COLVARS_OK;
  }
  int get_dE_dlambda(cvm::real *value) override
  {
    ++derivative_calls;
    *value = 0.0;
    return COLVARS_OK;
  }
  int send_alch_lambda() override { ++writeback_calls; return COLVARS_OK; }
  void reset_buffers()
  {
    energy = 0.0;
    for (auto &force : *modify_atom_applied_forces()) force.reset();
  }
  void distance(double value)
  {
    auto &positions = *modify_atom_positions();
    positions[0] = cvm::rvector(0.0, 0.0, 0.0);
    positions[1] = cvm::rvector(value, 0.0, 0.0);
  }
};

// Legal protected-member pointers allow inspection without casting a base
// object to a derived object or changing the production API for testing.
struct trial_probe : colvar {
  static colvarvalue const &extended(colvar const &cv) { return cv.*&trial_probe::x_ext; }
  static colvarvalue const &extended_velocity(colvar const &cv) { return cv.*&trial_probe::v_ext; }
  static bool restarted(colvar const &cv) { return cv.*&trial_probe::after_restart; }
  static double spring_energy(colvar const &cv) { return cv.*&trial_probe::potential_energy; }
  static double mass(colvar const &cv) { return cv.*&trial_probe::ext_mass; }
};

std::string configuration(bool extended, std::string const &extra = "")
{
  return "colvar {\n name d\n width 0.1\n" +
    std::string(extended ? " extendedLagrangian on\n extendedFluctuation 0.1\n extendedTimeConstant 100\n extendedLangevinDamping 1\n" : "") +
    extra +
    " distance {\n group1 { atomNumbers 1 }\n group2 { atomNumbers 2 }\n }\n}\n";
}

void configure(trial_proxy &proxy, std::string const &config)
{
  require(proxy.cvmodule->read_config_string(config) == COLVARS_OK, "configuration failed");
  proxy.distance(1.0);
  require(proxy.cvmodule->variables_active()->empty(), "initial active list unexpectedly populated");
}

void run_trial(trial_proxy &proxy)
{
  proxy.reset_buffers();
  require(proxy.cvmodule->calc_trial() == COLVARS_OK, "trial failed");
  require(!proxy.cvmodule->in_trial(), "trial flag leaked");
}

void test_plain_initial_trial()
{
  trial_proxy proxy;
  configure(proxy, configuration(false) +
    "harmonic {\n colvars d\n centers 0.5\n forceConstant 2\n}\n");
  auto *cv = proxy.cvmodule->colvar_by_name("d");
  double const value = cv->value().real_value;
  require(proxy.cvmodule->prepare_initial_trial_state() == COLVARS_OK, "prepare failed");
  run_trial(proxy);
  // Harmonic k is scaled by width^-2: k = 200.
  near(proxy.energy, 25.0, "pre-first-calc bias energy was not evaluated");
  near((*proxy.get_atom_applied_forces())[0].x, 100.0, "pre-first-calc bias force was not evaluated");
  near(cv->value().real_value, value, "trial changed reported value");
  require(proxy.cvmodule->variables_active()->empty(), "trial committed active list");
}

void test_extended_initial_trial(bool accept)
{
  trial_proxy proxy;
  configure(proxy, configuration(true));
  auto *cv = proxy.cvmodule->colvar_by_name("d");
  require(trial_probe::extended(*cv).type() == colvarvalue::type_notset, "fresh auxiliary coordinate initialized early");
  double const reported = cv->value().real_value;
  require(proxy.cvmodule->prepare_initial_trial_state() == COLVARS_OK, "prepare failed");
  double const k = cv->force_constant();
  for (double distance : {1.0, 1.1, 0.9, 1.1}) {
    proxy.distance(distance);
    run_trial(proxy);
    near(proxy.energy, 0.5 * k * (distance - 1.0) * (distance - 1.0), "initial trial recentered auxiliary coordinate");
    near((*proxy.get_atom_applied_forces())[0].x, k * (distance - 1.0), "initial spring force incorrect");
    require(trial_probe::extended(*cv).type() == colvarvalue::type_notset, "trial failed to restore uninitialized auxiliary state");
    near(cv->value().real_value, reported, "trial changed reported value");
  }
  require(proxy.random_calls == 0, "initial trial consumed random numbers");

  proxy.throw_energy = true;
  try {
    run_trial(proxy);
    throw std::logic_error("injected trial exception was not propagated");
  } catch (std::runtime_error const &error) {
    require(std::string(error.what()) == "injected trial failure", "wrong trial exception");
  }
  proxy.throw_energy = false;
  require(!proxy.cvmodule->in_trial(), "exception leaked trial flag");
  require(trial_probe::extended(*cv).type() == colvarvalue::type_notset, "exception leaked auxiliary overlay");

  double const distance = accept ? 1.1 : 1.0;
  proxy.distance(distance);
  proxy.reset_buffers();
  require(proxy.cvmodule->calc() == COLVARS_OK, "first real calculation failed");
  near(cv->value().real_value, 1.0, "first real calculation recentered accepted initial state");
  near((*proxy.get_atom_applied_forces())[0].x, k * (distance - 1.0), "first accepted/rejected real spring force incorrect");
  double const spring = 0.5 * k * (distance - 1.0) * (distance - 1.0);
  near(trial_probe::spring_energy(*cv), spring, "first accepted/rejected real spring energy incorrect");
  // The unchanged real integrator reports kinetic energy after its half kick.
  double const half_velocity = 0.5 * k * (distance - 1.0) / trial_probe::mass(*cv);
  near(proxy.energy, spring + 0.5 * trial_probe::mass(*cv) * half_velocity * half_velocity,
       "first real energy does not match the unchanged half-kick convention");
  require(proxy.random_calls == 1, "first real step did not integrate exactly once");
  proxy.cvmodule->clear_initial_trial_state();
}

void test_restart_and_boundaries()
{
  {
    trial_proxy proxy;
    configure(proxy, configuration(true));
    auto *cv = proxy.cvmodule->colvar_by_name("d");
    require(cv->set_state_params("name d\n x 1.0\n extended_x 0.95\n extended_v 0.002\n") == COLVARS_OK, "restart load failed");
    std::string const state = cv->get_state_params();
    require(proxy.cvmodule->prepare_initial_trial_state() == COLVARS_OK, "restart prepare failed");
    run_trial(proxy);
    require(cv->get_state_params() == state, "restart trial changed loaded state");
    require(trial_probe::restarted(*cv), "trial cleared restart validation");
    near(trial_probe::extended_velocity(*cv).real_value, 0.002, "restart velocity reset");
    proxy.reset_buffers();
    require(proxy.cvmodule->calc() == COLVARS_OK, "restart real calculation failed");
    near(cv->value().real_value, 0.95, "restart auxiliary coordinate recentered");
    near(cv->velocity().real_value, 0.002, "restart initial velocity reset");
    proxy.cvmodule->clear_initial_trial_state();
  }
  for (bool lower : {true, false}) {
    trial_proxy proxy;
    configure(proxy, configuration(true,
      " lowerBoundary 0.9\n upperBoundary 1.1\n reflectingLowerBoundary on\n reflectingUpperBoundary on\n"));
    double const center = lower ? 0.9 : 1.1;
    double const distance = lower ? 0.8 : 1.2;
    proxy.distance(distance);
    require(proxy.cvmodule->prepare_initial_trial_state() == COLVARS_OK, "reflecting prepare failed");
    auto *cv = proxy.cvmodule->colvar_by_name("d");
    run_trial(proxy);
    near(proxy.energy, 0.5 * cv->force_constant() * (distance-center)*(distance-center), "reflecting seed not clamped");
    proxy.reset_buffers();
    require(proxy.cvmodule->calc() == COLVARS_OK, "reflecting real calculation failed");
    near(cv->value().real_value, center, "real initialization did not use clamped seed");
    proxy.cvmodule->clear_initial_trial_state();
  }
}

void test_external_and_preparation_exception()
{
  trial_proxy proxy;
  require(proxy.cvmodule->read_config_string(
    "colvar {\n name lambda\n alchLambda {}\n extendedLagrangian on\n extendedMass 1\n extendedLangevinDamping 1\n}\n") == COLVARS_OK, "external configuration failed");
  auto *cv = proxy.cvmodule->colvar_by_name("lambda");
  proxy.indirect_lambda_biasing_force = 5.0;
  proxy.throw_lambda = true;
  try {
    proxy.cvmodule->prepare_initial_trial_state();
    throw std::logic_error("injected preparation exception was not propagated");
  } catch (std::runtime_error const &error) {
    require(std::string(error.what()) == "injected preparation failure", "wrong preparation exception");
  }
  require(!proxy.cvmodule->in_trial(), "preparation exception leaked trial flag");
  require(trial_probe::extended(*cv).type() == colvarvalue::type_notset, "preparation exception initialized state");
  proxy.throw_lambda = false;
  proxy.lambda = 0.42;
  require(proxy.cvmodule->prepare_initial_trial_state() == COLVARS_OK, "external prepare recovery failed");
  run_trial(proxy);
  require(proxy.derivative_calls == 0 && proxy.writeback_calls == 0 && proxy.random_calls == 0,
          "external preparation/trial read derivatives, wrote back, or consumed randomness");
  near(proxy.indirect_lambda_biasing_force, 5.0, "external trial cleared indirect force");
  require(trial_probe::extended(*cv).type() == colvarvalue::type_notset, "external trial leaked auxiliary state");
  proxy.lambda = 0.63;
  proxy.reset_buffers();
  require(proxy.cvmodule->calc() == COLVARS_OK, "external first real calculation failed");
  near(cv->value().real_value, 0.42, "external seed did not use initial Context parameter");
  proxy.cvmodule->clear_initial_trial_state();
}
}

int main()
{
  try {
    test_plain_initial_trial();
    test_extended_initial_trial(true);
    test_extended_initial_trial(false);
    test_restart_and_boundaries();
    test_external_and_preparation_exception();
    std::cout << "PASS: initial trial state, accept/reject, restart, boundaries, external parameter and exceptions\n";
  } catch (std::exception const &error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
