"""Analytical and integration regressions for the real OpenMM Force interface.

The oracle uses explicit mass-weighted COMs and a Gaussian/harmonic potential,
independently of the Colvars implementation. Units are nm, ps, Da, kJ/mol.
"""

import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import openmm as mm
from openmm import unit
import pytest

from openmm_colvars import ColvarsForce


MASSES = np.array([12.0, 11.0, 24.0, 36.0, 7.0, 48.0])
GROUPS = (np.array([0, 3]), np.array([2, 5]))
WIDTH, CENTER, FORCE_CONSTANT = 0.2, 0.35, 8.0
DT = 0.001


def distance_config(bias="harmonic", trajectory_frequency=0):
    common = f"""
units gromacs
colvarsTrajFrequency {trajectory_frequency}
colvarsRestartFrequency 0
colvar {{
  name d
  width {WIDTH}
  distance {{
    group1 {{ atomNumbers 1 4 }}
    group2 {{ atomNumbers 3 6 }}
  }}
}}
"""
    if bias == "harmonic":
        return common + f"""
harmonic {{
  name restraint
  colvars d
  centers {CENTER}
  forceConstant {FORCE_CONSTANT}
}}
"""
    if bias == "metadynamics":
        return common + """
metadynamics {
  name meta
  colvars d
  hillWeight 2.0
  gaussianSigmas 0.1
  newHillFrequency 1
  stepZeroData on
  useGrids off
  writeHillsTrajectory on
}
"""
    raise ValueError(bias)


def geometry(displacement):
    """Keep each group whole, vary the intergroup COM displacement only."""
    result = np.array([[0.18, 0.3, 0.4], [1.3, 0.8, 0.6],
                       [0.17, 0.3, 0.4], [0.24, 0.3, 0.4],
                       [0.6, 0.5, 0.9], [0.26, 0.3, 0.4]])
    com1 = np.average(result[GROUPS[0]], axis=0, weights=MASSES[GROUPS[0]])
    com2 = np.average(result[GROUPS[1]], axis=0, weights=MASSES[GROUPS[1]])
    result[GROUPS[1]] += com1 + np.array(displacement) - com2
    return result


def displacement_and_distance(positions, box=None):
    centers = [np.average(positions[g], axis=0, weights=MASSES[g]) for g in GROUPS]
    delta = centers[1] - centers[0]
    if box is not None:
        # Exhaustive nearby images: independent of the proxy's PBC code.
        delta = min((delta - np.array(image) @ np.asarray(box)
                     for image in itertools.product(range(-2, 3), repeat=3)),
                    key=lambda v: float(v @ v))
    return delta, float(np.linalg.norm(delta))


def harmonic_oracle(positions, box=None):
    delta, distance = displacement_and_distance(positions, box)
    spring = FORCE_CONSTANT / WIDTH**2
    energy = 0.5 * spring * (distance - CENTER)**2
    group_force = spring * (distance - CENTER) * delta / distance
    forces = np.zeros_like(positions)
    for group, sign in zip(GROUPS, (1.0, -1.0)):
        forces[group] = sign * (MASSES[group] / MASSES[group].sum())[:, None] * group_force
    return distance, energy, forces


def make_context(engine, tmp_path, config=None, box=None, group=1, from_file=False, extra_forces=()):
    system = mm.System()
    for mass in MASSES:
        system.addParticle(float(mass))
    for other in extra_forces:
        system.addForce(other)
    if box is not None:
        system.setDefaultPeriodicBoxVectors(*(mm.Vec3(*v) for v in box))
    text = distance_config() if config is None else config
    prefix = str(tmp_path / "run")
    if from_file:
        file = tmp_path / "colvars.in"
        file.write_text(text)
        force = ColvarsForce(str(file), prefix)
    else:
        force = ColvarsForce.fromConfigString(text, prefix)
    force.setUseGpu(engine.use_gpu)
    force.setUsesPeriodicBoundaryConditions(box is not None)
    force.setForceGroup(group)
    system.addForce(force)
    integrator = mm.VerletIntegrator(DT)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName(engine.platform), engine.properties)
    context.setVelocities(np.zeros((len(MASSES), 3)))
    return context, integrator, force


def state_values(context, groups=-1):
    state = context.getState(getEnergy=True, getForces=True, groups=groups)
    return (state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
            np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)))


def assert_analytical(context, force, engine, positions, box=None):
    expected_cv, expected_energy, expected_forces = harmonic_oracle(positions, box)
    energy, forces = state_values(context)
    cv = force.getColvarValue(context, "d")
    cv_tol, energy_tol, force_tol = engine.tolerances
    np.testing.assert_allclose(cv, expected_cv, rtol=0, atol=cv_tol)
    np.testing.assert_allclose(energy, expected_energy, rtol=0, atol=energy_tol)
    np.testing.assert_allclose(forces, expected_forces, rtol=0, atol=force_tol)
    np.testing.assert_allclose(forces.sum(axis=0), 0.0, rtol=0, atol=force_tol)
    np.testing.assert_array_equal(forces[[1, 4]], 0.0)
    print(json.dumps({"platform": engine.platform, "gpu": engine.use_gpu,
                      "precision": engine.precision, "backend": force.getBackend(context),
                      "cv_nm": cv, "expected_cv_nm": expected_cv,
                      "energy_kj_mol": energy, "expected_energy_kj_mol": expected_energy,
                      "max_force_error_kj_mol_nm": float(np.max(np.abs(forces - expected_forces))),
                      "forces_kj_mol_nm": forces.tolist(),
                      "absolute_tolerances": [cv_tol, energy_tol, force_tol]}))
    return energy, forces


@pytest.mark.parametrize("displacement,box", [
    ([0.7, 0.2, -0.1], None),
    ([1.8, 0.1, 0.0], [[2., 0, 0], [0, 2.2, 0], [0, 0, 2.4]]),
    ([1.8, 2.0, 0.1], [[2., 0, 0], [0.4, 2.2, 0], [-0.2, 0.3, 2.4]]),
], ids=("nonperiodic", "orthorhombic-boundary", "triclinic-boundary"))
def test_distance_harmonic_analytical(engine, tmp_path, displacement, box):
    context, integrator, force = make_context(engine, tmp_path, box=box, from_file=True)
    positions = geometry(displacement)
    context.setPositions(positions)
    assert_analytical(context, force, engine, positions, box)
    backend = force.getBackend(context).lower()
    if engine.platform == "CUDA" and engine.use_gpu:
        assert "colvarsgpu" in backend.replace(" ", "")
    if engine.platform == "CUDA" and not engine.use_gpu:
        assert "colvarscpu" in backend.replace(" ", "")
    del context, integrator


def test_bias_force_is_integrated(engine, tmp_path):
    context, integrator, force = make_context(engine, tmp_path)
    positions = geometry([0.7, 0.2, -0.1])
    context.setPositions(positions)
    _, _, expected_forces = harmonic_oracle(positions)
    integrator.step(1)
    observed = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    expected = positions + DT**2 * expected_forces / MASSES[:, None]
    tolerance = 2.0e-7 if engine.platform == "CUDA" and engine.precision == "single" else 2.0e-10
    np.testing.assert_allclose(observed, expected, rtol=0, atol=tolerance)
    del context, integrator


def test_repeated_queries_refresh_positions_and_box(engine, tmp_path):
    box = [[2., 0, 0], [0, 2.2, 0], [0, 0, 2.4]]
    context, integrator, force = make_context(engine, tmp_path, box=box)
    positions = geometry([1.8, 0.1, 0.0])
    context.setPositions(positions)
    energy, forces = assert_analytical(context, force, engine, positions, box)
    updates = force.getBiasUpdateCount(context)
    for _ in range(3):
        energy_only = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        force_only = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        np.testing.assert_allclose(energy_only, energy, rtol=0, atol=engine.tolerances[1])
        np.testing.assert_allclose(force_only, forces, rtol=0, atol=engine.tolerances[2])
    assert force.getBiasUpdateCount(context) == updates

    positions[GROUPS[1], 1] += 0.12
    context.setPositions(positions)
    assert_analytical(context, force, engine, positions, box)
    assert force.getBiasUpdateCount(context) == updates

    box[0][0] = 2.1
    context.setPeriodicBoxVectors(*(mm.Vec3(*v) for v in box))
    assert_analytical(context, force, engine, positions, box)
    assert force.getBiasUpdateCount(context) == updates
    del context, integrator


def test_tiny_coordinate_change_is_not_cached_away(engine, tmp_path):
    if engine.platform == "CUDA" and engine.precision == "single":
        pytest.skip("Sub-float-ulp displacement is not representable in CUDA single precision")
    context, integrator, force = make_context(engine, tmp_path)
    positions = geometry([0.7, 0.2, -0.1])
    context.setPositions(positions)
    first, _ = state_values(context)
    positions[GROUPS[1], 0] += 2.0e-7
    context.setPositions(positions)
    second, _ = assert_analytical(context, force, engine, positions)
    assert second != first
    del context, integrator


def test_force_group_queries_do_not_advance_bias(engine, tmp_path):
    context, integrator, force = make_context(engine, tmp_path, group=5)
    positions = geometry([0.7, 0.2, -0.1])
    context.setPositions(positions)
    excluded_energy, excluded_forces = state_values(context, groups=1 << 0)
    assert excluded_energy == 0.0
    np.testing.assert_array_equal(excluded_forces, 0.0)
    assert force.getBiasUpdateCount(context) == 0
    energy, forces = state_values(context, groups=1 << 5)
    _, expected_energy, expected_forces = harmonic_oracle(positions)
    np.testing.assert_allclose(energy, expected_energy, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(forces, expected_forces, rtol=0, atol=engine.tolerances[2])
    assert force.getBiasUpdateCount(context) == 1
    del context, integrator


def test_metadynamics_same_step_is_read_only(engine, tmp_path):
    context, integrator, force = make_context(engine, tmp_path, config=distance_config("metadynamics"))
    positions = geometry([0.6, 0.0, 0.0])
    context.setPositions(positions)
    initial_energy, initial_forces = state_values(context)
    np.testing.assert_allclose(initial_energy, 2.0, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(initial_forces, 0.0, rtol=0, atol=engine.tolerances[2])
    count = force.getBiasUpdateCount(context)
    for _ in range(4):
        repeated_energy, _ = state_values(context)
        np.testing.assert_allclose(repeated_energy, initial_energy, rtol=0, atol=engine.tolerances[1])
    assert force.getBiasUpdateCount(context) == count

    positions[GROUPS[1], 0] += 0.05
    context.setPositions(positions)
    energy, forces = state_values(context)
    expected_energy = 2.0 * math.exp(-0.5 * (0.05 / 0.1)**2)
    np.testing.assert_allclose(energy, expected_energy, rtol=0, atol=engine.tolerances[1])
    expected = np.zeros_like(positions)
    # Positive displacement from a hill pushes the groups farther apart.
    for group, sign in zip(GROUPS, (-1.0, 1.0)):
        expected[group, 0] = sign * expected_energy * 0.05 / 0.1**2 * MASSES[group] / MASSES[group].sum()
    np.testing.assert_allclose(forces, expected, rtol=0, atol=engine.tolerances[2])
    assert force.getBiasUpdateCount(context) == count
    # Advancing the OpenMM step permits exactly one new state update.
    integrator.step(1)
    state_values(context)
    assert force.getBiasUpdateCount(context) == count + 1
    del context, integrator


@pytest.mark.parametrize("bad_index", (0, 7, -1))
def test_invalid_atom_index_rejected(engine, tmp_path, bad_index):
    config = distance_config().replace("atomNumbers 1 4", f"atomNumbers {bad_index}")
    with pytest.raises(Exception, match="(?i)atom|index|particle"):
        make_context(engine, tmp_path, config=config)


def test_unit_mismatch_rejected(engine, tmp_path):
    with pytest.raises(Exception, match="(?i)unit|gromacs"):
        make_context(engine, tmp_path, config=distance_config().replace("units gromacs", "units real"))


def test_bias_adds_to_physical_forces(engine, tmp_path):
    physical = mm.CustomExternalForce("0.5*k*(x*x+y*y+z*z)")
    physical.addGlobalParameter("k", 3.0)
    for index in range(len(MASSES)):
        physical.addParticle(index, [])
    physical.setForceGroup(0)
    context, integrator, force = make_context(engine, tmp_path, extra_forces=(physical,))
    positions = geometry([0.7, 0.2, -0.1])
    context.setPositions(positions)
    _, bias_energy, bias_forces = harmonic_oracle(positions)
    physical_energy = 1.5 * float(np.sum(positions**2))
    physical_forces = -3.0 * positions
    energy, forces = state_values(context)
    np.testing.assert_allclose(energy, bias_energy + physical_energy, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(forces, bias_forces + physical_forces, rtol=0, atol=engine.tolerances[2])
    physical_energy_only, physical_forces_only = state_values(context, groups=1)
    np.testing.assert_allclose(physical_energy_only, physical_energy, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(physical_forces_only, physical_forces, rtol=0, atol=engine.tolerances[2])
    assert force.getBiasUpdateCount(context) == 1
    integrator.step(1)
    observed = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    expected = positions + DT**2 * (bias_forces + physical_forces) / MASSES[:, None]
    tolerance = 2.0e-7 if engine.platform == "CUDA" and engine.precision == "single" else 2.0e-10
    np.testing.assert_allclose(observed, expected, rtol=0, atol=tolerance)
    del context, integrator


def test_two_force_instances_have_independent_state(engine, tmp_path):
    system = mm.System()
    for mass in MASSES:
        system.addParticle(float(mass))
    forces = []
    for index in range(2):
        force = ColvarsForce.fromConfigString(distance_config(), str(tmp_path / f"force{index}"))
        force.setUseGpu(engine.use_gpu)
        force.setUsesPeriodicBoundaryConditions(False)
        force.setForceGroup(index + 1)
        system.addForce(force)
        forces.append(force)
    integrator = mm.VerletIntegrator(DT)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName(engine.platform), engine.properties)
    positions = geometry([0.7, 0.2, -0.1])
    context.setPositions(positions)
    energy, atomic_forces = state_values(context)
    _, one_energy, one_forces = harmonic_oracle(positions)
    np.testing.assert_allclose(energy, 2.0 * one_energy, rtol=0, atol=2 * engine.tolerances[1])
    np.testing.assert_allclose(atomic_forces, 2.0 * one_forces, rtol=0, atol=2 * engine.tolerances[2])
    assert [f.getBiasUpdateCount(context) for f in forces] == [1, 1]
    context.setStepCount(1)
    state_values(context, groups=1 << 1)
    assert [f.getBiasUpdateCount(context) for f in forces] == [2, 1]
    state_values(context, groups=1 << 2)
    assert [f.getBiasUpdateCount(context) for f in forces] == [2, 2]
    del context, integrator


def test_moving_restraint_does_not_repeat_work_or_center_update(engine, tmp_path):
    config = distance_config().replace("forceConstant 8.0", """forceConstant 8.0
  targetCenters 0.55
  targetNumSteps 10
  outputAccumulatedWork on""")
    context, integrator, force = make_context(engine, tmp_path, config=config)
    positions = geometry([0.7, 0.0, 0.0])
    context.setPositions(positions)
    state_values(context)
    context.setStepCount(1)
    energy, _ = state_values(context)
    state = force.getBiasState(context)
    center = float(re.search(r"\bcenters\s+([-+0-9.eE]+)", state).group(1))
    work = float(re.search(r"\baccumulatedWork\s+([-+0-9.eE]+)", state).group(1))
    np.testing.assert_allclose(center, 0.37, rtol=0, atol=1e-12)
    np.testing.assert_allclose(energy, 0.5 * FORCE_CONSTANT / WIDTH**2 * (0.7 - center)**2,
                               rtol=0, atol=engine.tolerances[1])
    count = force.getBiasUpdateCount(context)
    positions[GROUPS[1], 0] += 0.03
    context.setPositions(positions)
    energy, _ = state_values(context)
    for _ in range(3):
        state_values(context)
    repeated = force.getBiasState(context)
    assert float(re.search(r"\bcenters\s+([-+0-9.eE]+)", repeated).group(1)) == center
    assert float(re.search(r"\baccumulatedWork\s+([-+0-9.eE]+)", repeated).group(1)) == work
    np.testing.assert_allclose(energy, 0.5 * FORCE_CONSTANT / WIDTH**2 * (0.73 - center)**2,
                               rtol=0, atol=engine.tolerances[1])
    assert force.getBiasUpdateCount(context) == count
    del context, integrator


def test_well_tempered_hill_height_and_repeat_safety(engine, tmp_path):
    config = distance_config("metadynamics").replace("useGrids off", """useGrids off
  wellTempered on
  biasTemperature 1200""")
    context, integrator, force = make_context(engine, tmp_path, config=config)
    positions = geometry([0.6, 0.0, 0.0])
    context.setPositions(positions)
    first, _ = state_values(context)
    np.testing.assert_allclose(first, 2.0, rtol=0, atol=engine.tolerances[1])
    context.setStepCount(1)
    second, _ = state_values(context)
    second_hill = 2.0 * math.exp(-2.0 / (0.00831446261815324 * 1200.0))
    np.testing.assert_allclose(second, 2.0 + second_hill, rtol=0, atol=engine.tolerances[1])
    count = force.getBiasUpdateCount(context)
    positions[GROUPS[1], 0] += 0.05
    context.setPositions(positions)
    expected = (2.0 + second_hill) * math.exp(-0.5 * (0.05 / 0.1)**2)
    for _ in range(3):
        energy, _ = state_values(context)
        np.testing.assert_allclose(energy, expected, rtol=0, atol=engine.tolerances[1])
    assert force.getBiasUpdateCount(context) == count
    del context, integrator


def test_rewinding_step_counter_is_rejected(engine, tmp_path):
    context, integrator, force = make_context(engine, tmp_path)
    context.setPositions(geometry([0.7, 0.0, 0.0]))
    state_values(context)
    integrator.step(2)
    state_values(context)
    context.setStepCount(0)
    with pytest.raises(Exception, match="(?i)backwards|restart|step count"):
        state_values(context)
    del context, integrator


def test_overlapping_context_lifetimes_are_safe(engine, tmp_path):
    # Colvars historically has static feature descriptors. A subprocess turns
    # a lifetime regression into a clear failure instead of aborting pytest.
    code = f"""
import gc
import sys
from pathlib import Path
sys.path.insert(0, {str(Path(__file__).resolve().parent)!r})
from conftest import Engine
from test_openmm_colvars import make_context, geometry, state_values, distance_config
engine = Engine({engine.platform!r}, {engine.use_gpu!r}, {engine.precision!r})
root = Path({str(tmp_path)!r})
for name in ('first', 'second', 'bad'):
    (root/name).mkdir()
first = make_context(engine, root/'first')
second = make_context(engine, root/'second')
for context, integrator, force in (first, second):
    context.setPositions(geometry([0.7, 0.0, 0.0]))
    state_values(context)
del context, integrator, force
try:
    make_context(engine, root/'bad', config=distance_config().replace('atomNumbers 1 4', 'atomNumbers 9999'))
except Exception:
    pass
else:
    raise AssertionError('invalid constructor unexpectedly succeeded')
del first
gc.collect()
context, integrator, force = second
integrator.step(1)
state_values(context)
assert force.getBiasUpdateCount(context) == 2
del context, integrator, force, second
gc.collect()
print('CONTEXT_LIFETIMES_OK')
"""
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout[-6000:] + result.stderr[-6000:]
    assert "CONTEXT_LIFETIMES_OK" in result.stdout


def test_coordnum_pairlist_refreshes_after_coordinate_and_box_edits(engine, tmp_path):
    config = """
units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {
  name contact
  width 1.0
  coordNum {
    group1 { atomNumbers 1 }
    group2 { atomNumbers 3 }
    cutoff 0.4
    expNumer 6
    expDenom 12
    tolerance 0.01
    pairListFrequency 100
  }
}
harmonic {
  name contacts
  colvars contact
  centers 0.5
  forceConstant 3.0
}
"""
    box = np.diag([12., 12., 12.])
    context, integrator, force = make_context(engine, tmp_path, config=config, box=box)
    positions = np.zeros((len(MASSES), 3))

    def check(x, step, side):
        context.setStepCount(step)
        context.setPeriodicBoxVectors(mm.Vec3(side, 0, 0), mm.Vec3(0, side, 0), mm.Vec3(0, 0, side))
        positions[2, 0] = x
        context.setPositions(positions)
        energy, forces = state_values(context)
        signed_distance = x - side * round(x / side)
        r = abs(signed_distance)
        switching = 1.0 / (1.0 + (r / 0.4)**6)
        expected_cv = max(0.0, (switching - 0.01) / 0.99)
        expected_energy = 1.5 * (expected_cv - 0.5)**2
        expected_forces = np.zeros_like(positions)
        if expected_cv > 0.0:
            derivative = -(6.0 * r**5 / 0.4**6) * switching**2 / 0.99
            expected_forces[0, 0] = 3.0 * (expected_cv - 0.5) * derivative * np.sign(signed_distance)
            expected_forces[2, 0] = -expected_forces[0, 0]
        np.testing.assert_allclose(force.getColvarValue(context, "contact"), expected_cv,
                                   rtol=0, atol=engine.tolerances[0])
        np.testing.assert_allclose(energy, expected_energy, rtol=0, atol=engine.tolerances[1])
        np.testing.assert_allclose(forces, expected_forces, rtol=0, atol=engine.tolerances[2])

    check(5.0, 0, 12.0)  # Build an empty neighbor list.
    check(0.3, 1, 12.0)  # Enter range before pairListFrequency steps elapse.
    check(5.0, 1, 12.0)
    check(0.3, 1, 12.0)  # Also invalidate within the same OpenMM step.
    check(5.0, 1, 12.0)
    check(5.0, 1, 5.3)   # A box edit alone turns a far pair into a near pair.
    assert force.getBiasUpdateCount(context) == 2
    del context, integrator
