"""Regression tests for OpenMM reference coordinates and extended dynamics."""

import math
import os
from pathlib import Path
import tempfile
import unittest

import openmm as mm
from openmm import unit
from openmm_colvars import ColvarsForce


PLATFORMS = os.environ.get("OPENMM_COLVARS_TEST_PLATFORMS", "Reference").split()
TETRAHEDRON = [(0.1, 0.1, 0.1), (0.1, -0.1, -0.1), (-0.1, 0.1, -0.1), (-0.1, -0.1, 0.1)]
RMSD_CONFIG = """
colvarsTrajFrequency 1
colvar {{
  name r
  width 1.0
  rmsd {{
    atoms {{ atomNumbers 1 2 3 4 }}
    {reference}
  }}
}}
harmonic {{
  colvars r
  centers 0.0
  forceConstant 1.0
}}
"""
EABF_CONFIG = """
colvarsTrajFrequency 1
colvar {{
  name d
  width 0.1
  lowerBoundary 0.5
  upperBoundary 1.5
  outputVelocity on
  extendedLagrangian on
  extendedTemp 300
  extendedFluctuation 0.05
  extendedTimeConstant 200
  extendedLangevinDamping {damping}
  distance {{
    group1 {{ atomNumbers 1 }}
    group2 {{ atomNumbers 2 }}
  }}
}}
abf {{
  colvars d
  applyBias off
}}
"""


def make_context(config, prefix, platform_name, positions, integrator, tabf_readback=False, alchemy=False):
    system = mm.System()
    for _ in positions:
        system.addParticle(12.0)
    force = ColvarsForce.fromConfigString(config, str(prefix))
    force.setRandomSeed(17)
    force.setTabfEnabled(tabf_readback)
    force.setAlchemyEnabled(alchemy)
    system.addForce(force)
    platform = mm.Platform.getPlatformByName(platform_name)
    properties = {"Precision": "double"} if platform_name == "CUDA" else {}
    context = mm.Context(system, integrator, platform, properties)
    context.setPositions([mm.Vec3(*xyz) for xyz in positions] * unit.nanometer)
    return system, context


def run(config, prefix, platform_name, positions, dt_fs=1, steps=0,
        initial_dt_fs=None, change_after_steps=0, probe_before_change=False, tabf_readback=False):
    # Evaluate forces at every integration step.
    start_dt_fs = dt_fs if initial_dt_fs is None else initial_dt_fs
    integrator = mm.CustomIntegrator(start_dt_fs * unit.femtoseconds)
    integrator.addPerDofVariable("observedForce", 0)
    integrator.addComputePerDof("observedForce", "f")
    if tabf_readback:
        # Translate both atoms for the first two steps so the TABF cache is
        # refreshed. Keep the last geometry unchanged for the final readback.
        integrator.addGlobalVariable("phase", 0)
        integrator.addComputePerDof("x", "x+0.001*step(1.5-phase)")
        integrator.addComputeGlobal("phase", "phase+1")
    else:
        integrator.addComputePerDof("x", "x")  # Invalidate OpenMM's force cache.
    system, context = make_context(config, prefix, platform_name, positions, integrator, tabf_readback)
    if initial_dt_fs is not None:
        if change_after_steps:
            integrator.step(change_after_steps)
        if probe_before_change:
            context.getState(energy=True, forces=True)
        integrator.setStepSize(dt_fs * unit.femtoseconds)
    if steps > change_after_steps:
        integrator.step(steps - change_after_steps)
    state = context.getState(energy=True, forces=True)
    del context  # Close the Colvars trajectory before reading it.
    lines = prefix.with_suffix(".colvars.traj").read_text().splitlines()
    columns = lines[0].split()[1:]
    return state, [dict(zip(columns, map(float, line.split()))) for line in lines[1:]]


class TestUnitsAndTimestep(unittest.TestCase):
    def test_rmsd_reference_coordinate_units(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xyz = root / "ref.xyz"
            xyz.write_text("4\nReference in angstrom\n" + "".join(
                f"C {10*x:.3f} {10*y:.3f} {10*z:.3f}\n" for x, y, z in TETRAHEDRON
            ))
            pdb = root / "ref.pdb"
            pdb.write_text("".join(
                f"ATOM  {i:5d}  C{i:<2d} UNK A   1    {10*x:8.3f}{10*y:8.3f}{10*z:8.3f}"
                f"{1.0:6.2f}{0.0:6.2f}           C\n"
                for i, (x, y, z) in enumerate(TETRAHEDRON, 1)
            ) + "END\n")
            inline = "refPositions " + " ".join(f"({x}, {y}, {z})" for x, y, z in TETRAHEDRON)
            references = {"inline": inline, "xyz": f"refPositionsFile {xyz}", "pdb": f"refPositionsFile {pdb}"}
            for platform in PLATFORMS:
                for mode, reference in references.items():
                    with self.subTest(platform=platform, mode=mode):
                        config = RMSD_CONFIG.format(reference=reference)
                        state, trajectory = run(config, root / f"{platform}-{mode}", platform, TETRAHEDRON)
                        self.assertAlmostEqual(trajectory[0]["r"], 0.0, delta=1e-6)
                        self.assertAlmostEqual(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),
                                               0.0, delta=1e-10)
                        for force in state.getForces().value_in_unit(unit.kilojoules_per_mole / unit.nanometer):
                            self.assertLess(math.sqrt(sum(component**2 for component in force)), 1e-6)

    def test_eabf_noise_uses_integrator_step_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            positions = [(0, 0, 0), (1, 0, 0)]
            for platform in PLATFORMS:
                with self.subTest(platform=platform):
                    first_kick = {}
                    for dt_fs, damping in [(1, 1), (2, 1), (2, 0)]:
                        config = EABF_CONFIG.format(damping=damping)
                        prefix = root / f"{platform}-{dt_fs}fs-{damping}"
                        _, trajectory = run(config, prefix, platform, positions, dt_fs=dt_fs, steps=1)
                        # At step zero the extended coordinate is at the physical value and at rest.
                        # Its step-one velocity is therefore entirely the first random kick.
                        self.assertEqual(trajectory[0]["vr_d"], 0.0)
                        first_kick[dt_fs, damping] = trajectory[1]["vr_d"]
                    self.assertNotEqual(first_kick[1, 1], 0.0)
                    self.assertEqual(first_kick[2, 0], 0.0)
                    # The same seed gives the same Gaussian deviate in both runs.
                    expected_ratio = math.sqrt(math.expm1(-0.004) / math.expm1(-0.002))
                    self.assertAlmostEqual(first_kick[2, 1] / first_kick[1, 1], expected_ratio, delta=1e-8)

    def test_eabf_step_size_change_after_force_probe(self):
        config = EABF_CONFIG.format(damping=1)
        positions = [(0, 0, 0), (1, 0, 0)]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for platform in PLATFORMS:
                for boundary in (0, 2):
                    with self.subTest(platform=platform, change_after_steps=boundary):
                        options = dict(dt_fs=2, steps=4, initial_dt_fs=1, change_after_steps=boundary)
                        prefix = root / f"{platform}-{boundary}"
                        _, expected = run(config, prefix, platform, positions, **options)
                        _, observed = run(config, prefix.with_name(prefix.name + "-probe"), platform, positions,
                                          probe_before_change=True, **options)
                        self.assertEqual([r["step"] for r in expected], list(range(5)))
                        self.assertEqual([r["step"] for r in observed], list(range(5)))
                        for result, baseline in zip(observed, expected):
                            self.assertAlmostEqual(result["r_d"], baseline["r_d"], delta=1e-11)
                            self.assertAlmostEqual(result["vr_d"], baseline["vr_d"], delta=1e-11)

                        if boundary:
                            # Query the cached energy while the extended spring is stretched.
                            options["steps"] = boundary
                            state, _ = run(config, root / f"{platform}-energy", platform, positions, **options)
                            probe_state, _ = run(config, root / f"{platform}-probe-energy", platform, positions,
                                                 probe_before_change=True, **options)
                            self.assertAlmostEqual(
                                probe_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),
                                state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),
                                delta=1e-12,
                            )

    def test_tabf_post_step_readback_keeps_completed_step(self):
        config = EABF_CONFIG.format(damping=1)
        positions = [(0, 0, 0), (1, 0, 0)]
        with tempfile.TemporaryDirectory() as tmp:
            for platform in PLATFORMS:
                with self.subTest(platform=platform):
                    options = dict(steps=3, initial_dt_fs=1, change_after_steps=3, tabf_readback=True)
                    prefix = Path(tmp) / platform
                    baseline, trajectory = run(config, prefix, platform, positions, dt_fs=1, **options)
                    changed, recorded = run(config, prefix.with_name(platform + "-changed"),
                                            platform, positions, dt_fs=2, **options)
                    # TABF permits reusing the preceding step for a readback of the same geometry.
                    self.assertEqual([r["step"] for r in trajectory], [0, 1, 2])
                    self.assertEqual(recorded, trajectory)
                    self.assertEqual(changed.getPotentialEnergy(), baseline.getPotentialEnergy())

    def test_adaptive_steps_are_rejected_for_extended_dynamics(self):
        config = EABF_CONFIG.format(damping=1)
        positions = [(0, 0, 0), (1, 0, 0)]
        with tempfile.TemporaryDirectory() as tmp:
            for platform in PLATFORMS:
                for integrator_name in ("VariableVerlet", "VariableLangevin", "Compound"):
                    with self.subTest(platform=platform, integrator=integrator_name):
                        if integrator_name == "VariableVerlet":
                            integrator = mm.VariableVerletIntegrator(0.001)
                        elif integrator_name == "VariableLangevin":
                            integrator = mm.VariableLangevinIntegrator(300, 1, 0.001)
                        else:
                            integrator = mm.CompoundIntegrator()
                            integrator.addIntegrator(mm.VerletIntegrator(0.001))
                            integrator.addIntegrator(mm.VariableVerletIntegrator(0.001))
                        prefix = Path(tmp) / f"{platform}-{integrator_name}"
                        system, context = make_context(config, prefix, platform, positions, integrator)
                        try:
                            if integrator_name == "Compound":
                                context.getState(forces=True)
                                integrator.setCurrentIntegrator(1)
                            with self.assertRaisesRegex(mm.OpenMMException, "adaptive integrators are not supported"):
                                context.getState(forces=True)
                        finally:
                            del context

    def test_step_size_change_with_changed_positions_is_rejected(self):
        positions = [(0, 0, 0), (1, 0, 0)]
        geometry = EABF_CONFIG.format(damping=1)
        external = geometry + """
colvar {
  name l
  extendedLagrangian on
  extendedMass 100
  extendedTemp 300
  extendedLangevinDamping 0
  alchLambda { }
}
"""
        with tempfile.TemporaryDirectory() as tmp:
            for platform in PLATFORMS:
                for alchemy, config, error in (
                    (False, geometry, "reusing the previous force evaluation"),
                    (True, external, "extended external coordinate"),
                ):
                    with self.subTest(platform=platform, alchemy=alchemy):
                        integrator = mm.VerletIntegrator(0.001)
                        prefix = Path(tmp) / f"{platform}-{alchemy}"
                        system, context = make_context(config, prefix, platform, positions, integrator, alchemy=alchemy)
                        try:
                            if alchemy:
                                context.setParameter("lambda", 0.5)
                            context.getState(forces=True)
                            integrator.setStepSize(0.002)
                            context.setPositions([mm.Vec3(0, 0, 0), mm.Vec3(1.01, 0, 0)])
                            with self.assertRaisesRegex(mm.OpenMMException, error):
                                context.getState(forces=True)
                        finally:
                            del context

    def test_failed_timestep_replay_does_not_return_stale_cache(self):
        config = (EABF_CONFIG.format(damping=1)
                  .replace("width 0.1", "width 0.00002")
                  .replace("lowerBoundary 0.5", "lowerBoundary 0.99995\n  reflectingLowerBoundary on")
                  .replace("upperBoundary 1.5", "upperBoundary 1.00005\n  reflectingUpperBoundary on"))
        with tempfile.TemporaryDirectory() as tmp:
            for platform in PLATFORMS:
                with self.subTest(platform=platform):
                    integrator = mm.VerletIntegrator(0.001)
                    system, context = make_context(config, Path(tmp) / platform, platform,
                                                   [(0, 0, 0), (1, 0, 0)], integrator)
                    try:
                        context.getState(energy=True, forces=True)
                        integrator.setStepSize(0.004)
                        with self.assertRaisesRegex(mm.OpenMMException, "still outside boundaries"):
                            context.getState(energy=True, forces=True)
                        integrator.setStepSize(0.001)
                        with self.assertRaisesRegex(mm.OpenMMException, "still outside boundaries"):
                            context.getState(energy=True, forces=True)
                    finally:
                        del context


if __name__ == "__main__":
    unittest.main()
