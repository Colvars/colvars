"""Regression tests for the force snapshot supplied to ordinary OpenMM ABF."""

import os
from pathlib import Path
import tempfile
import unittest

import openmm as mm
from openmm import unit
from openmm_colvars import ColvarsForce


CONFIG = """
colvarsTrajFrequency 1
colvar {
  name d
  width 0.1
  lowerBoundary 0.5
  upperBoundary 1.5
  outputTotalForce on
  distance {
    group1 { atomNumbers 1 }
    group2 { atomNumbers 2 }
  }
}
abf {
  colvars d
  applyBias off
  outputFreq 1
}
"""


FORCE_UNIT = unit.kilojoules_per_mole / unit.nanometer
PLATFORMS = os.environ.get("OPENMM_COLVARS_TEST_PLATFORMS", "Reference").split()


def table(path):
    return [
        [float(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def make_system(mode):
    system = mm.System()
    for mass in (12.0, 12.0, 12.0 if mode == "no_sites" else 0.0, 12.0):
        system.addParticle(mass)
    if mode != "no_sites":
        w = 1.0 / 3.0 if mode == "rounding" else 0.25
        system.setVirtualSite(2, mm.TwoParticleAverageSite(0, 1, w, 1.0 - w))
    if mode in ("nonbonded", "no_sites", "pme"):
        force = mm.NonbondedForce()
        if mode == "pme":
            system.setDefaultPeriodicBoxVectors(mm.Vec3(3, 0, 0), mm.Vec3(0, 3, 0), mm.Vec3(0, 0, 3))
            force.setNonbondedMethod(mm.NonbondedForce.PME)
            force.setCutoffDistance(0.75)
            force.setReciprocalSpaceForceGroup(3)
        for i in range(4):
            charged = 1 if mode == "no_sites" else 2
            force.addParticle(1.0 if i == charged else -1.0 if i == 3 else 0.0, 1.0, 0.0)
    else:
        force = mm.CustomExternalForce("-fx*x")
        force.addPerParticleParameter("fx")
        if mode == "nested":
            site = system.addParticle(0.0)
            system.setVirtualSite(site, mm.TwoParticleAverageSite(0, 2, 0.5, 0.5))
            force.addParticle(site, [-8.0])
        elif mode == "direct":
            force.addParticle(0, [1.0])
            force.addParticle(1, [3.0])
        else:
            force.addParticle(2, [1.0e6 if mode == "rounding" else 4.0])
    system.addForce(force)
    return system


class TestABFVirtualSites(unittest.TestCase):
    def run_system(self, platform_name, mode, energy_queries=False, partial_force_query=False,
                   apply_bias=False):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "abf"
            system = make_system(mode)
            config = CONFIG
            if energy_queries:
                config += """
                histogram {
                  name hist
                  colvars d
                  outputFreq 1
                }
                """
            if apply_bias:
                config = config.replace("applyBias off", "applyBias on\n  fullSamples 1")
            colvars = ColvarsForce.fromConfigString(config, str(prefix))
            colvars.setForceGroup(2)
            # Remove the kBT-dependent Jacobian force in this mechanical test.
            colvars.setTargetTemperature(0.0)
            system.addForce(colvars)
            integrator = mm.CustomIntegrator(0.001)
            integrator.setKineticEnergyExpression("0.5*m*v*v")
            platform = mm.Platform.getPlatformByName(platform_name)
            properties = {"Threads": "2"} if platform_name == "CPU" else {}
            if mode == "rounding" and platform_name == "CUDA":
                properties["Precision"] = "mixed"
            context = mm.Context(system, integrator, platform, properties)
            xyz = [0.25, 1.25, 1.0, 2.0]
            if mode == "nested":
                xyz.append(0.625)
            context.setPositions([mm.Vec3(x, 0.0, 0.0) for x in xyz])
            context.computeVirtualSites()
            # Keep positions fixed so the expected force and sampled bin are constant.
            # Step zero initializes Colvars and is not counted by default.
            for step in range(4):
                context.setStepCount(step)
                if energy_queries and step in (1, 2):
                    context.getState(energy=True)
                if energy_queries and step == 2:
                    continue
                if partial_force_query and step == 1:
                    context.getState(forces=True, groups={0, 2} if mode == "pme" else {2})
                # An explicit mask of all active groups is equivalent to the default mask.
                if partial_force_query and step == 3:
                    state = context.getState(forces=True, groups={0, 2, 3} if mode == "pme" else {0, 2})
                else:
                    state = context.getState(forces=True)
            forces = state.getForces().value_in_unit(FORCE_UNIT)
            del context, integrator, system, colvars

            traj = Path(f"{prefix}.colvars.traj")
            labels = traj.read_text().splitlines()[0].lstrip("#").split()
            ft = labels.index("ft_d")
            force_rows = [(int(row[0]), row[ft]) for row in table(traj)]
            counts = dict(table(Path(f"{prefix}.count")))
            grads = dict(table(Path(f"{prefix}.grad")))
            occupied = [(n, grads[x]) for x, n in counts.items() if n > 0]
            self.assertEqual(len(occupied), 1)
            if energy_queries:
                histogram = table(Path(f"{prefix}.hist.dat"))
                self.assertEqual(sum(row[1] for row in histogram), 2.0)
            return forces, force_rows, occupied[0]

    def test_equivalent_forces(self):
        for platform in PLATFORMS:
            for mode, parent_x in (("direct", (1.0, 3.0)), ("site", (1.0, 3.0)),
                                   ("nested", (-5.0, -3.0))):
                with self.subTest(platform=platform, mode=mode):
                    forces, rows, (count, grad) = self.run_system(platform, mode)
                    self.assertEqual(tuple(v[0] for v in forces[:2]), parent_x)
                    self.assertEqual(rows, [(step, 1.0) for step in range(4)])
                    self.assertEqual((count, grad), (3.0, -1.0))

    def test_native_nonbonded_forces(self):
        for platform in PLATFORMS:
            for mode in ("nonbonded", "no_sites"):
                with self.subTest(platform=platform, mode=mode):
                    forces, rows, (count, grad) = self.run_system(platform, mode)
                    expected = (forces[1][0] - forces[0][0]) / 2
                    self.assertGreater(expected, 30.0)
                    for _, projected in rows:
                        self.assertAlmostEqual(projected, expected, delta=1e-6)
                    self.assertEqual(count, 3)
                    self.assertAlmostEqual(grad, -expected, delta=1e-6)

    def test_snapshot_matches_native_virtual_site_rounding(self):
        for platform in PLATFORMS:
            with self.subTest(platform=platform):
                forces, rows, (count, grad) = self.run_system(platform, "rounding")
                expected = (forces[1][0] - forces[0][0]) / 2
                self.assertGreater(expected, 1.0e5)
                for _, projected in rows:
                    self.assertAlmostEqual(projected, expected, delta=1e-7)
                self.assertEqual(count, 3)
                self.assertAlmostEqual(grad, -expected, delta=1e-7)

    def test_partial_force_query_does_not_replace_real_sample(self):
        for platform in PLATFORMS:
            with self.subTest(platform=platform):
                _, rows, (count, grad) = self.run_system(platform, "site", partial_force_query=True)
                self.assertEqual(rows, [(step, 1.0) for step in range(4)])
                self.assertEqual((count, grad), (3.0, -1.0))

    def test_pme_reciprocal_group_is_required_for_abf_sampling(self):
        for platform in PLATFORMS:
            with self.subTest(platform=platform):
                forces, rows, (count, grad) = self.run_system(platform, "pme", partial_force_query=True)
                expected = (forces[1][0] - forces[0][0]) / 2
                self.assertGreater(expected, 10.0)
                self.assertEqual([step for step, _ in rows], list(range(4)))
                for _, projected in rows:
                    self.assertAlmostEqual(projected, expected, delta=1e-6)
                self.assertEqual(count, 3)
                self.assertAlmostEqual(grad, -expected, delta=1e-6)

    def test_bias_forces_are_added_to_the_live_buffer(self):
        for platform in PLATFORMS:
            with self.subTest(platform=platform):
                forces, rows, (count, grad) = self.run_system(platform, "site", apply_bias=True)
                self.assertEqual(rows, [(step, 1.0) for step in range(4)])
                self.assertEqual((count, grad), (3.0, -1.0))
                self.assertEqual(tuple(v[0] for v in forces[:2]), (2.0, 2.0))

    def test_energy_queries_do_not_record_atomic_force_samples(self):
        for platform in PLATFORMS:
            with self.subTest(platform=platform):
                _, rows, (count, grad) = self.run_system(platform, "site", energy_queries=True)
                self.assertEqual((count, grad), (2.0, -1.0))
                self.assertEqual(rows, [(step, 1.0) for step in (0, 1, 3)])


if __name__ == "__main__":
    unittest.main()
