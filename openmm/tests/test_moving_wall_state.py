"""Moving-restraint snapshots restore the instantaneous bias and its schedule."""

import re

import numpy as np
import openmm as mm
import pytest

from openmm_colvars import ColvarsForce
from test_openmm_colvars import state_values


def configuration(kind="upper", staged=False, work=True):
    bias_type = "harmonicWalls"
    walls = "upperWalls 0.5\n  targetUpperWalls 0.7"
    if kind == "lower":
        walls = "lowerWalls 0.5\n  targetLowerWalls 0.7"
    elif kind == "both":
        walls = "lowerWalls 0.3\n  targetLowerWalls 0.2\n  upperWalls 0.7\n  targetUpperWalls 0.8"
    elif kind == "upper-and-k":
        walls += "\n  targetForceConstant 2.0"
    elif kind == "spring-k":
        bias_type = "harmonic"
        walls = "centers 0.5\n  targetForceConstant 2.0"
    steps = "targetNumSteps 4\n  targetNumStages 3" if staged else "targetNumSteps 10"
    work_option = "outputAccumulatedWork on" if work and not staged else ""
    return f"""units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {{
  name d
  width 0.1
  distance {{
    group1 {{ atomNumbers 1 }}
    group2 {{ atomNumbers 2 }}
  }}
}}
{bias_type} {{
  name moving
  colvars d
  {walls}
  {steps}
  forceConstant 1
  {work_option}
}}
"""


def make_context(engine, output, config, distance=1.0, state_file=None):
    force = ColvarsForce.fromConfigString(config, str(output))
    force.setUseGpu(engine.use_gpu)
    force.setUsesPeriodicBoundaryConditions(False)
    if state_file is not None:
        force.setRestartInputPrefix(str(state_file))
    system = mm.System()
    system.addParticle(12)
    system.addParticle(24)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName(engine.platform), engine.properties)
    context.setPositions([[0, 0, 0], [distance, 0, 0]])
    return context, integrator, force


def schedule_state(force, context):
    text = force.getBiasState(context)
    result = {}
    for key in ("firstStep", "stage", "dAdLambda", "accumulatedWork",
                "lowerWalls", "upperWalls", "forceConstant"):
        values = re.findall(r"^\s*" + key + r"\s+([^\n]+)$", text, re.MULTILINE)
        if values:
            assert len(values) == 1, f"Duplicate {key} in state: {text}"
            result[key] = tuple(float(x) for x in values[0].split())
    return result


def assert_equal(actual, expected, engine):
    np.testing.assert_allclose(actual[0], expected[0], rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(actual[1], expected[1], rtol=0, atol=engine.tolerances[2])


def assert_schedule_equal(actual, expected):
    assert actual.keys() == expected.keys()
    for key in expected:
        np.testing.assert_allclose(actual[key], expected[key], rtol=0, atol=1e-11,
                                   err_msg=f"Restart changed {key}")


@pytest.mark.parametrize("kind", ("upper", "lower", "both", "upper-and-k", "spring-k"))
@pytest.mark.parametrize("save_step", (5, 10))
def test_continuous_moving_restraint_restart(engine, tmp_path, kind, save_step):
    config = configuration(kind)
    distance = 0.1 if kind == "lower" else 1.0
    original, integrator, force = make_context(engine, tmp_path / "original", config, distance)
    for step in range(save_step + 1):
        original.setStepCount(step)
        expected = state_values(original)
    schedule = schedule_state(force, original)
    assert "accumulatedWork" in schedule
    assert schedule["firstStep"] == (0.0,)
    if kind == "upper" and save_step == 10:
        np.testing.assert_allclose(expected[0], 4.5, rtol=0, atol=engine.tolerances[1])
        np.testing.assert_allclose(expected[1], [[30, 0, 0], [-30, 0, 0]],
                                   rtol=0, atol=engine.tolerances[2])
    state_file = tmp_path / "saved.colvars.state"
    state_file.write_text(force.getBiasState(original))
    restored, restored_integrator, restored_force = make_context(
        engine, tmp_path / "restored", config, distance, state_file)
    restored.setStepCount(37)
    for _ in range(3):
        assert_equal(state_values(restored), expected, engine)
        assert_schedule_equal(schedule_state(restored_force, restored), schedule)
    assert restored_force.getBiasUpdateCount(restored) == 0

    # Changed geometry uses the frozen restored potential without adding work.
    changed = [[0, 0, 0], [distance + 0.015, 0, 0]]
    original.setPositions(changed)
    restored.setPositions(changed)
    assert_equal(state_values(restored), state_values(original), engine)
    assert_schedule_equal(schedule_state(restored_force, restored), schedule)
    for context in (original, restored):
        context.setPositions([[0, 0, 0], [distance, 0, 0]])
    del context

    for step in range(save_step + 1, 13):
        original.setStepCount(step)
        restored.setStepCount(37 + step - save_step)
        assert_equal(state_values(restored), state_values(original), engine)
        assert_schedule_equal(schedule_state(restored_force, restored), schedule_state(force, original))
    del restored, restored_integrator, original, integrator


@pytest.mark.parametrize("kind", ("upper", "both", "spring-k"))
def test_staged_restraint_restores_partial_stage_sum(engine, tmp_path, kind):
    config = configuration(kind, staged=True)
    original, integrator, force = make_context(engine, tmp_path / "original", config)
    for step in range(6):
        original.setStepCount(step)
        expected = state_values(original)
    schedule = schedule_state(force, original)
    assert schedule["stage"] == (1.0,)
    assert schedule["dAdLambda"][0] != 0.0
    state_file = tmp_path / "saved.colvars.state"
    state_file.write_text(force.getBiasState(original))
    restored, restored_integrator, restored_force = make_context(
        engine, tmp_path / "restored", config, state_file=state_file)
    for _ in range(3):
        assert_equal(state_values(restored), expected, engine)
        assert_schedule_equal(schedule_state(restored_force, restored), schedule)
    assert restored_force.getBiasUpdateCount(restored) == 0
    for step in range(6, 14):
        original.setStepCount(step)
        restored.setStepCount(step - 5)
        assert_equal(state_values(restored), state_values(original), engine)
        assert_schedule_equal(schedule_state(restored_force, restored), schedule_state(force, original))
    del restored, restored_integrator, original, integrator


def test_nonzero_wall_start_step_is_preserved(engine, tmp_path):
    config = configuration()
    original, integrator, force = make_context(engine, tmp_path / "original", config)
    for step in range(6):
        original.setStepCount(step)
        expected = state_values(original)
    snapshot = force.getBiasState(original)
    # The same physical schedule can have been introduced seven steps later.
    # Shift both saved clocks together, retaining the instantaneous wall/work.
    snapshot = re.sub(r"^(\s*step\s+)(\d+)",
                      lambda m: m[1] + str(int(m[2]) + 7), snapshot, flags=re.MULTILINE)
    snapshot = re.sub(r"^(\s*firstStep\s+)0\b", r"\g<1>7", snapshot, flags=re.MULTILINE)
    state_file = tmp_path / "shifted.colvars.state"
    state_file.write_text(snapshot)
    restored, restored_integrator, restored_force = make_context(
        engine, tmp_path / "restored", config, state_file=state_file)
    assert_equal(state_values(restored), expected, engine)
    assert schedule_state(restored_force, restored)["firstStep"] == (7.0,)
    for step in range(6, 13):
        original.setStepCount(step)
        restored.setStepCount(step - 5)
        assert_equal(state_values(restored), state_values(original), engine)
        np.testing.assert_allclose(schedule_state(restored_force, restored)["accumulatedWork"],
                                   schedule_state(force, original)["accumulatedWork"], rtol=0, atol=1e-11)
    del restored, restored_integrator, original, integrator


@pytest.mark.parametrize("known_start", (False, True))
def test_legacy_wall_state_requires_reconstructible_schedule(engine, tmp_path, known_start):
    config = configuration(work=False)
    original, integrator, force = make_context(engine, tmp_path / "original", config)
    for step in range(11):
        original.setStepCount(step)
        expected = state_values(original)
    snapshot = force.getBiasState(original)
    snapshot = re.sub(r"^\s*upperWalls[^\n]*\n", "", snapshot, flags=re.MULTILINE)
    if not known_start:
        snapshot = re.sub(r"^\s*firstStep[^\n]*\n", "", snapshot, flags=re.MULTILINE)
    state_file = tmp_path / "legacy.colvars.state"
    state_file.write_text(snapshot)
    if known_start:
        restored, restored_integrator, restored_force = make_context(
            engine, tmp_path / "restored", config, state_file=state_file)
        assert_equal(state_values(restored), expected, engine)
        restored.setStepCount(1)
        assert_equal(state_values(restored), expected, engine)
        del restored, restored_integrator
    else:
        with pytest.raises(Exception, match="(?i)firstStep|legacy.*schedule"):
            make_context(engine, tmp_path / "rejected", config, state_file=state_file)
    del original, integrator


@pytest.mark.parametrize("save_step", (4, 5))
def test_legacy_staged_state_needs_partial_sum_inside_stage(engine, tmp_path, save_step):
    config = configuration("spring-k", staged=True)
    original, integrator, force = make_context(engine, tmp_path / "original", config)
    for step in range(save_step + 1):
        original.setStepCount(step)
        expected = state_values(original)
    snapshot = re.sub(r"^\s*dAdLambda[^\n]*\n", "", force.getBiasState(original), flags=re.MULTILINE)
    state_file = tmp_path / "legacy.colvars.state"
    state_file.write_text(snapshot)
    if save_step == 4:
        restored, restored_integrator, restored_force = make_context(
            engine, tmp_path / "restored", config, state_file=state_file)
        assert_equal(state_values(restored), expected, engine)
        assert schedule_state(restored_force, restored)["dAdLambda"] == (0.0,)
        del restored, restored_integrator
    else:
        with pytest.raises(Exception, match="dAdLambda"):
            make_context(engine, tmp_path / "rejected", config, state_file=state_file)
    del original, integrator


def test_single_moving_wall_does_no_work_inside_the_wall(engine, tmp_path):
    context, integrator, force = make_context(engine, tmp_path / "inside", configuration(), distance=0.2)
    for step in range(12):
        context.setStepCount(step)
        energy, forces = state_values(context)
        np.testing.assert_allclose(energy, 0.0, rtol=0, atol=engine.tolerances[1])
        np.testing.assert_allclose(forces, 0.0, rtol=0, atol=engine.tolerances[2])
        assert schedule_state(force, context)["accumulatedWork"] == (0.0,)
    del context, integrator


@pytest.mark.parametrize("kind,distance", [
    ("upper", 1.0), ("lower", 0.1), ("both", 1.0), ("both", 0.1),
    ("upper-and-k", 1.0), ("static-lower", 0.1), ("upper", 0.2),
])
def test_staged_wall_derivative_matches_potential_finite_difference(engine, tmp_path, kind, distance):
    config = configuration("upper" if kind == "static-lower" else kind, staged=True)
    if kind == "static-lower":
        config = config.replace("upperWalls 0.5", "lowerWalls 0.3\n  upperWalls 0.5")
    context, integrator, force = make_context(engine, tmp_path / "derivative", config, distance)

    def potential(lam):
        lower, upper = None, None
        if kind == "lower":
            lower = 0.5 + 0.2 * lam
        elif kind == "both":
            lower, upper = 0.3 - 0.1 * lam, 0.7 + 0.1 * lam
        else:
            upper = 0.5 + 0.2 * lam
            if kind == "static-lower":
                lower = 0.3
        displacement = 0.0
        if lower is not None and distance < lower:
            displacement = distance - lower
        if upper is not None and distance > upper:
            displacement = distance - upper
        spring = 1.0 + lam if kind == "upper-and-k" else 1.0
        return 0.5 * spring / 0.1**2 * displacement**2

    # After step 1, the partial stage sum contains exactly one derivative at
    # lambda=0. After the boundary at step4 it is reset; step5 contains exactly
    # one derivative at lambda=1/3. This tests the accumulated observable through
    # the public state API, without calling the implementation's derivative.
    for step in range(6):
        context.setStepCount(step)
        state_values(context)
        if step in (1, 5):
            lam = 0.0 if step == 1 else 1.0 / 3.0
            eps = 1e-5
            numerical_derivative = (potential(lam + eps) - potential(lam - eps)) / (2 * eps)
            actual = schedule_state(force, context)["dAdLambda"][0]
            np.testing.assert_allclose(actual, numerical_derivative, rtol=0,
                                       atol=max(1e-7, engine.tolerances[2]))
    del context, integrator
