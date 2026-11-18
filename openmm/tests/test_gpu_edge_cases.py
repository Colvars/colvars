"""Analytical regressions also exercised by the CUDA GPU and CPU-fallback modes."""

import numpy as np
import openmm as mm
from openmm import unit
import pytest

from openmm_colvars import ColvarsForce


@pytest.mark.parametrize("dummy_first", [False, True])
def test_distance_to_dummy_group(engine, tmp_path, dummy_first):
    # Only particle 2 is selected, so the dummy group must not accidentally
    # create an atom slot or alter the force on the unrelated first particle.
    dummy = np.array([0.1, 0.2, 0.3])
    groups = ["atomNumbers 2", "dummyAtom (0.1, 0.2, 0.3)"]
    if dummy_first:
        groups.reverse()
    config = f"""
units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {{
  name d
  width 0.5
  distance {{
    group1 {{ {groups[0]} }}
    group2 {{ {groups[1]} }}
  }}
}}
harmonic {{
  colvars d
  centers 0.4
  forceConstant 10.0
}}
"""
    system = mm.System()
    system.addParticle(7.0)
    system.addParticle(12.0)
    force = ColvarsForce.fromConfigString(config, str(tmp_path / "dummy"))
    force.setUseGpu(engine.use_gpu)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName(engine.platform), engine.properties)
    for selected_position in ([0.7, 0.8, 0.9], [0.8, -0.1, 0.7]):
        context.setPositions([[9.0, 8.0, 7.0], selected_position])
        displacement = np.asarray(selected_position) - dummy
        distance = float(np.linalg.norm(displacement))
        spring = 10.0 / 0.5**2
        expected_energy = 0.5 * spring * (distance - 0.4)**2
        expected_force = -spring * (distance - 0.4) * displacement / distance
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        observed = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        np.testing.assert_allclose(force.getColvarValue(context, "d"), distance,
                                   atol=engine.tolerances[0], rtol=0)
        np.testing.assert_allclose(energy, expected_energy,
                                   atol=engine.tolerances[1], rtol=0)
        np.testing.assert_allclose(observed[1], expected_force,
                                   atol=engine.tolerances[2], rtol=0)
        np.testing.assert_array_equal(observed[0], 0.0)
    assert force.getBiasUpdateCount(context) == 1
    del context, integrator


def test_only_dummy_groups_request_no_atoms(engine, tmp_path):
    config = """
units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {
  name d
  width 0.5
  distance {
    group1 { dummyAtom (0.1, 0.2, 0.3) }
    group2 { dummyAtom (0.8, -0.1, 0.7) }
  }
}
harmonic {
  colvars d
  centers 0.4
  forceConstant 10.0
}
"""
    system = mm.System()
    system.addParticle(12.0)
    force = ColvarsForce.fromConfigString(config, str(tmp_path / "all_dummy"))
    force.setUseGpu(engine.use_gpu)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName(engine.platform), engine.properties)
    context.setPositions([[9.0, 8.0, 7.0]])
    distance = float(np.linalg.norm([0.7, -0.3, 0.4]))
    expected_energy = 0.5 * (10.0 / 0.5**2) * (distance - 0.4)**2
    state = context.getState(getEnergy=True, getForces=True)
    np.testing.assert_allclose(force.getColvarValue(context, "d"), distance,
                               atol=engine.tolerances[0], rtol=0)
    np.testing.assert_allclose(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                               expected_energy, atol=engine.tolerances[1], rtol=0)
    np.testing.assert_array_equal(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer),
                                  0.0)
    if engine.platform == "CUDA":
        assert "selected atoms=0" in force.getBackend(context)
    del context, integrator
