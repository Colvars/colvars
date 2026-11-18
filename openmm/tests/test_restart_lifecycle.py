"""A loaded adaptive step must be evaluated without repeating its update."""

import re

import numpy as np
import openmm as mm
import pytest

from openmm_colvars import ColvarsForce
from test_openmm_colvars import state_values


CONFIG = """
units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {
  name d
  width 0.1
  lowerBoundary 0
  upperBoundary 2
  distance {
    group1 { atomNumbers 1 }
    group2 { atomNumbers 2 }
  }
}
metadynamics {
  name meta
  colvars d
  hillWeight 2
  gaussianSigmas 0.1
  newHillFrequency 1
  stepZeroData on
  useGrids on
  gridsUpdateFrequency 1
}
"""


def context_with_restart(engine, output, state_file=None, load_source="prefix"):
    config = CONFIG
    if state_file is not None and load_source == "configuration":
        config += f"\ndefaultInputStateFile {state_file}\n"
    force = ColvarsForce.fromConfigString(config, str(output))
    force.setUseGpu(engine.use_gpu)
    force.setUsesPeriodicBoundaryConditions(False)
    if state_file is not None and load_source == "prefix":
        force.setRestartInputPrefix(str(state_file))
    system = mm.System()
    system.addParticle(12.0)
    system.addParticle(24.0)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName(engine.platform), engine.properties)
    context.setPositions([[0, 0, 0], [0.65, 0, 0]])
    return context, integrator, force


def assert_same_state(actual, expected, engine):
    np.testing.assert_allclose(actual[0], expected[0], rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(actual[1], expected[1], rtol=0, atol=engine.tolerances[2])


@pytest.mark.parametrize("load_source", ["prefix", "configuration"])
def test_restored_step_is_not_deposited_twice(engine, tmp_path, load_source):
    original, integrator, force = context_with_restart(engine, tmp_path / "original")
    for step in range(3):
        original.setStepCount(step)
        expected = state_values(original)
    np.testing.assert_allclose(expected[0], 6.0, rtol=0, atol=engine.tolerances[1])
    state_file = tmp_path / "saved.colvars.state"
    state_file.write_text(force.getBiasState(original))

    restored, restored_integrator, restored_force = context_with_restart(
        engine, tmp_path / "restored", state_file, load_source)
    # The caller may continue with an independently restored OpenMM step count.
    # The initial query still represents the saved Colvars step (2).
    restored.setStepCount(37)
    assert_same_state(state_values(restored), expected, engine)
    assert_same_state(state_values(restored), expected, engine)
    assert restored_force.getBiasUpdateCount(restored) == 0
    snapshot = restored_force.getBiasState(restored)
    assert int(re.search(r"configuration\s*\{\s*step\s+(\d+)", snapshot)[1]) == 2

    # This change exceeds the restart geometry-check threshold. Only the first
    # loaded geometry is checked; arbitrary later same-step queries are valid.
    changed_positions = [[0, 0, 0], [0.85, 0, 0]]
    original.setPositions(changed_positions)
    restored.setPositions(changed_positions)
    expected = state_values(original)
    assert_same_state(state_values(restored), expected, engine)
    assert restored_force.getBiasUpdateCount(restored) == 0

    original.setStepCount(3)
    restored.setStepCount(38)
    expected = state_values(original)
    assert_same_state(state_values(restored), expected, engine)
    assert_same_state(state_values(restored), expected, engine)
    assert restored_force.getBiasUpdateCount(restored) == 1
    snapshot = restored_force.getBiasState(restored)
    assert int(re.search(r"configuration\s*\{\s*step\s+(\d+)", snapshot)[1]) == 3
    del restored, restored_integrator, original, integrator
