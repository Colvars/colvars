"""Deterministic local tests of replica state/buffer replacement ordering.

These exercise filesystem communication in one process; they do not establish
the correctness or performance of a distributed filesystem or MPI transport.
"""

import re

import numpy as np
import openmm as mm
from openmm import unit
import pytest

from openmm_colvars import ColvarsForce


def make_replica(engine, directory, replica, *, grids):
    grid_options = ("gridsUpdateFrequency 1\n  keepHills off"
                    if grids else "")
    config = f"""
units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {{
  name d
  width 0.1
  lowerBoundary 0
  upperBoundary 4
  distance {{
    group1 {{ atomNumbers 1 }}
    group2 {{ atomNumbers 2 }}
  }}
}}
metadynamics {{
  name meta
  colvars d
  hillWeight 2
  gaussianSigmas 0.1
  newHillFrequency 1
  stepZeroData on
  useGrids {'on' if grids else 'off'}
  {grid_options}
  multipleReplicas on
  replicaID {replica}
  replicasRegistry {directory/'registry'}
  replicaUpdateFrequency 1
}}
"""
    force = ColvarsForce.fromConfigString(config, str(directory/replica))
    force.setUseGpu(engine.use_gpu)
    force.setUsesPeriodicBoundaryConditions(False)
    system = mm.System()
    system.addParticle(12.0)
    system.addParticle(12.0)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName(engine.platform), engine.properties)
    # An interior grid node: projected hills contribute their full height, and
    # keepHills-off snapshots contain no explicit or boundary hills.
    context.setPositions([[0, 0, 0], [2.05, 0, 0]])
    return context, integrator, force


def energy(context, step):
    context.setStepCount(step)
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


@pytest.mark.parametrize("grids", [False, True])
@pytest.mark.parametrize("legacy", [False, True])
def test_replica_snapshot_overlap_skips_old_buffer_but_accepts_new_hill(
        engine, tmp_path, monkeypatch, grids, legacy):
    monkeypatch.chdir(tmp_path)
    first = make_replica(engine, tmp_path, "A", grids=grids)
    for step in (0, 1):
        energy(first[0], step)
    snapshot = first[2].getBiasState(first[0])
    snapshot = snapshot[snapshot.index("metadynamics {"):]
    if grids:
        assert "hill {" not in snapshot
    else:
        assert snapshot.count("hill {") == 2
    assert "lastHillStep 1" in snapshot
    if legacy:
        snapshot = re.sub(r"lastHillStep [^\n]+\n", "", snapshot)

    # Reproduce the producer's real replacement window: its state file has
    # advanced, but its old incremental buffer has not yet been truncated.
    (tmp_path/"A.colvars.meta.A.state").write_text(snapshot)
    buffer = tmp_path/"A.colvars.meta.A.hills"
    assert buffer.read_text().count("hill {") == 2
    second = make_replica(engine, tmp_path, "B", grids=grids)
    np.testing.assert_allclose(energy(second[0], 0), 6.0, rtol=0, atol=engine.tolerances[1])

    # A genuinely newer record must still be accepted once, despite the old
    # records remaining in the same buffer. B adds its own hill at each step.
    np.testing.assert_allclose(energy(first[0], 2), 8.0, rtol=0, atol=engine.tolerances[1])
    assert buffer.read_text().count("hill {") == 3
    np.testing.assert_allclose(energy(second[0], 1), 10.0, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(energy(second[0], 1), 10.0, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(energy(second[0], 2), 12.0, rtol=0, atol=engine.tolerances[1])
    del first, second


@pytest.mark.parametrize("grids", [False, True])
def test_initial_empty_replica_snapshot_accepts_step_zero_hill(
        engine, tmp_path, monkeypatch, grids):
    monkeypatch.chdir(tmp_path)
    first = make_replica(engine, tmp_path, "A", grids=grids)
    initial_snapshot = (tmp_path/"A.colvars.meta.A.state").read_text()
    assert "lastHillStep -1" in initial_snapshot
    assert "hill {" not in initial_snapshot
    np.testing.assert_allclose(energy(first[0], 0), 2.0, rtol=0, atol=engine.tolerances[1])

    # The snapshot predates the first deposition even though both have module
    # step zero. Its explicit -1 watermark permits the new buffered step-zero
    # hill, which a plain h_it <= state_file_step test would discard.
    second = make_replica(engine, tmp_path, "B", grids=grids)
    np.testing.assert_allclose(energy(second[0], 0), 4.0, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(energy(first[0], 1), 6.0, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(energy(second[0], 1), 8.0, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(energy(second[0], 1), 8.0, rtol=0, atol=engine.tolerances[1])
    del first, second
