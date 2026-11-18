"""Reject configurations whose native Colvars step semantics OpenMM cannot honor."""

import pytest

from test_openmm_colvars import distance_config, make_context


@pytest.mark.parametrize("scope", ["colvar", "bias"])
def test_colvars_multiple_time_step_is_rejected(engine, tmp_path, scope):
    config = distance_config()
    if scope == "colvar":
        config = config.replace("name d\n", "name d\n  timeStepFactor 2\n")
    else:
        config = config.replace("name restraint\n", "name restraint\n  timeStepFactor 2\n")
    with pytest.raises(Exception, match="timeStepFactor|multiple.time"):
        make_context(engine, tmp_path, config=config)
