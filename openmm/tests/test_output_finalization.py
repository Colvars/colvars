"""End-of-run output must preserve buffered hills and the current potential."""

from pathlib import Path

import pytest

from test_metadynamics_state import configuration, make_context, evaluate, assert_same


@pytest.mark.parametrize("grids", [False, True])
@pytest.mark.parametrize("output_frequency", [0, 2])
def test_explicit_flush_saves_pending_output_without_updating_bias(
        engine, tmp_path, grids, output_frequency):
    config = configuration(grids=grids).replace(
        "newHillFrequency 1",
        f"newHillFrequency 1\n  writeHillsTrajectory on\n  outputFreq {output_frequency}")
    prefix = tmp_path / "output"
    original = make_context(engine, prefix, config)
    context, _, force = original
    for step in range(4):
        expected = evaluate(context, 0.65, step)
    snapshot = force.getBiasState(context)
    updates = force.getBiasUpdateCount(context)
    hills_file = Path(str(prefix) + ".colvars.meta.hills.traj")
    state_file = Path(str(prefix) + ".colvars.state")

    for _ in range(2):
        force.flushOutputs(context)
        steps = [int(line.split()[0]) for line in hills_file.read_text().splitlines()
                 if line.strip() and not line.startswith("#")]
        assert steps == [0, 1, 2, 3]
        assert state_file.read_text() == snapshot
        assert force.getBiasState(context) == snapshot
        assert force.getBiasUpdateCount(context) == updates
        assert_same(evaluate(context, 0.65), expected, engine)
    if grids:
        assert Path(str(prefix) + ".pmf").is_file()

    restored = make_context(engine, tmp_path / "restored", config, state_file)
    assert_same(evaluate(restored[0], 0.65, 3), expected, engine)
    # Flushing is repeatable and does not terminate or reset the simulation.
    assert_same(evaluate(restored[0], 0.65, 4), evaluate(context, 0.65, 4), engine)
    force.flushOutputs(context)
    steps = [int(line.split()[0]) for line in hills_file.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    assert steps == [0, 1, 2, 3, 4]


def test_flush_requires_an_evaluated_state(engine, tmp_path):
    original = make_context(engine, tmp_path / "empty", configuration(grids=False))
    context, _, force = original
    with pytest.raises(Exception, match="Evaluate energy or forces"):
        force.flushOutputs(context)
    evaluate(context, 0.65)
    force.flushOutputs(context)
