"""Check the shipped input against the two-particle analytical solution."""

import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit
from openmm_colvars import ColvarsForce


def test_minimal_example(engine, tmp_path):
    masses = np.array([12., 24.])
    initial = np.array([[0., 0., 0.], [0.65, 0., 0.]])
    expected_forces = np.array([[15., 0., 0.], [-15., 0., 0.]])
    dt = 0.001
    system = mm.System()
    for mass in masses:
        system.addParticle(mass)
    # SWIG passes Python paths as UTF-8, including on Windows legacy code pages.
    tmp_path = tmp_path / "偏置-Δ"
    tmp_path.mkdir()
    source = Path(__file__).resolve().parents[1] / "examples" / "colvars.in"
    config = tmp_path / "输入.in"
    config.write_bytes(source.read_bytes())
    force = ColvarsForce(str(config), str(tmp_path / "minimal"))
    force.setUseGpu(engine.use_gpu)
    force.setUsesPeriodicBoundaryConditions(False)
    force.setForceGroup(1)
    system.addForce(force)
    integrator = mm.VerletIntegrator(dt)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName(engine.platform),
                         engine.properties)
    context.setPositions(initial)
    context.setVelocities(np.zeros_like(initial))

    def check(quantity, actual, expected, tolerance):
        error = np.abs(np.asarray(actual) - expected)
        passed = bool(np.all(error <= tolerance))
        print(json.dumps(dict(case="two-particle harmonic", quantity=quantity,
                              platform=engine.platform, precision=engine.precision,
                              backend=force.getBackend(context), actual=np.asarray(actual).tolist(),
                              expected=np.asarray(expected).tolist(), absolute_error=error.tolist(),
                              tolerance=tolerance, status="PASS" if passed else "FAIL")))
        np.testing.assert_allclose(actual, expected, atol=tolerance, rtol=0)

    excluded = context.getState(getEnergy=True, getForces=True, groups={0})
    check("excluded energy (kJ/mol)", excluded.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), 0., 0.)
    check("excluded forces (kJ/mol/nm)", excluded.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer), np.zeros_like(initial), 0.)
    assert force.getBiasUpdateCount(context) == 0
    cv_tol, energy_tol, force_tol = engine.tolerances
    for _ in range(2):
        state = context.getState(getEnergy=True, getForces=True, groups={1})
        check("distance (nm)", force.getColvarValue(context, "distance"), 0.65, cv_tol)
        check("energy (kJ/mol)", state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
              1.125, energy_tol)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        check("forces (kJ/mol/nm)", forces, expected_forces, force_tol)
        check("net force (kJ/mol/nm)", forces.sum(axis=0), np.zeros(3), force_tol)
        assert force.getBiasUpdateCount(context) == 1

    integrator.step(1)
    final = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    check("positions after one step (nm)", final, initial + dt**2 * expected_forces / masses[:, None],
          2e-7 if engine.precision == "single" else 2e-10)
    context.getState(getEnergy=True)
    assert context.getStepCount() == 1 and force.getBiasUpdateCount(context) == 2
    force.flushOutputs(context)
    assert (tmp_path / "minimal.colvars.state").stat().st_size > 0
    assert (tmp_path / "minimal.colvars.traj").stat().st_size > 0
