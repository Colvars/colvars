"""OpenMM Monte Carlo volume trials and the Colvars callback lifecycle."""

import json

import numpy as np
import openmm as mm
from openmm import unit
import pytest

from test_moving_wall_state import configuration
from test_openmm_colvars import assert_analytical, distance_config, geometry, make_context


@pytest.mark.parametrize("config", [
    distance_config("metadynamics"),
    distance_config().replace("centers 0.35", "centers 0.35\n  targetCenters 0.6\n  targetNumSteps 10"),
    configuration("upper"),
    configuration("spring-k", staged=True),
], ids=["metadynamics", "moving-center", "moving-wall", "staged-force-constant"])
def test_adaptive_barostat_guard(engine, tmp_path, config):
    barostat = mm.MonteCarloBarostat(1.0, 300.0, 1)
    barostat.setName("pressure control")  # detection must survive user naming
    with pytest.raises(Exception, match="static restraints with a Monte Carlo barostat"):
        make_context(engine, tmp_path, config=config, box=np.eye(3) * 2, extra_forces=[barostat])


def test_static_restraint_accepted_and_rejected_volume_trials(engine, tmp_path):
    # Each integrator step attempts one volume move with no MD displacement.
    # High pressure favors contractions, so both outcomes occur frequently.
    barostat = mm.MonteCarloBarostat(100000.0, 300.0, 1)
    barostat.setRandomNumberSeed(137)
    trial_step = mm.CustomIntegrator(0.001)
    trial_step.addUpdateContextState()
    box = np.eye(3) * 2
    context, integrator, force = make_context(engine, tmp_path, box=box, extra_forces=[barostat],
                                              integrator=trial_step)
    context.setPositions(geometry([1.8, 0.1, 0.0]))
    before = context.getState(getPositions=True)
    accepted = rejected = 0
    for _ in range(64):
        integrator.step(1)
        after = context.getState(getPositions=True)
        positions = after.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        box = after.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        if after.getPeriodicBoxVolume() == before.getPeriodicBoxVolume():
            rejected += 1
            np.testing.assert_array_equal(positions, before.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        else:
            accepted += 1
        assert_analytical(context, force, engine, positions, box)
        before = after
        if accepted and rejected:
            break
    assert accepted and rejected
    print(json.dumps({"platform": engine.platform, "precision": engine.precision, "gpu": engine.use_gpu,
                      "accepted_volume_trials": accepted, "rejected_volume_trials": rejected}))
