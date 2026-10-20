"""Tests for ColvarsCollectiveVariable, a serializable force whose energy is the value
of a Colvars collective variable.

Run with: python -m unittest discover openmm/tests/python
"""

import copy
import math
import os
import pickle
import tempfile
import unittest

import numpy as np
import openmm
from openmm_colvars import ColvarsCollectiveVariable, ColvarsForce

PLATFORMS = ["Reference", "CPU"]
RADIANS_TO_DEGREES = 180 / math.pi


def make_positions(num_particles, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(scale=0.3, size=(num_particles, 3))


def make_system(num_particles, *forces):
    system = openmm.System()
    for _ in range(num_particles):
        system.addParticle(12.0)
    for group, force in enumerate(forces):
        force.setForceGroup(group)
        system.addForce(force)
    return system


def make_context(system, positions, platform="Reference"):
    context = openmm.Context(
        system,
        openmm.VerletIntegrator(0.001),
        openmm.Platform.getPlatformByName(platform),
    )
    context.setPositions(positions)
    return context


def energy_and_forces(context, group):
    state = context.getState(getEnergy=True, getForces=True, groups={group})
    energy = state.getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)
    forces = state.getForces(asNumpy=True).value_in_unit(
        openmm.unit.kilojoule_per_mole / openmm.unit.nanometer
    )
    return energy, forces


def colvar_config(name, component, extra=""):
    return f"colvar {{\n  name {name}\n{extra}  {component}\n}}\n"


def atoms(*indices):
    return " ".join(str(i + 1) for i in indices)


DISTANCE = (
    "distance {\n    group1 { atomNumbers "
    + atoms(0)
    + " }\n    group2 { atomNumbers "
    + atoms(1)
    + " }\n  }"
)
DIHEDRAL = (
    "dihedral {\n"
    + "".join(f"    group{k + 1} {{ atomNumbers {atoms(k)} }}\n" for k in range(4))
    + "  }"
)

# A customFunction converts the dihedral from degrees to radians
DIHEDRAL_IN_RADIANS = colvar_config(
    "phi",
    DIHEDRAL.replace("dihedral {\n", "dihedral {\n    name t\n", 1),
    "  customFunction t*0.017453292519943295\n",
)


def rotation_matrix(seed):
    rng = np.random.default_rng(seed)
    w, x, y, z = (q := rng.normal(size=4)) / np.linalg.norm(q)
    matrix = np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )
    return matrix, 2 * math.degrees(math.acos(min(1.0, abs(w))))


def orientation_angle(reference):
    ref_positions = "\n".join(
        f"      ({x:.10f}, {y:.10f}, {z:.10f})" for x, y, z in reference
    )
    component = (
        "orientationAngle {\n    atoms { atomNumbers "
        + atoms(*range(len(reference)))
        + " }\n    refPositions {\n"
        + ref_positions
        + "\n    }\n  }"
    )
    return colvar_config("omega", component)


def configured_cv():
    cv = ColvarsCollectiveVariable(colvar_config("phi", DIHEDRAL), "phi")
    cv.setForceGroup(7)
    cv.setName("phi_colvars")
    cv.setUsesPeriodicBoundaryConditions(True)
    return cv


def coupling_force(phi, distance, expression):
    # A periodic harmonic coupling of a dihedral in radians with a distance
    force = openmm.CustomCVForce(expression)
    force.addGlobalParameter("k", 100.0)
    force.addGlobalParameter("phi0", 2.5)
    force.addCollectiveVariable("phi", phi)
    force.addCollectiveVariable("d", distance)
    return force


class TestColvarsCollectiveVariable(unittest.TestCase):
    def assert_matches(self, context, rtol=1e-6):
        energy, forces = energy_and_forces(context, 0)
        expected_energy, expected_forces = energy_and_forces(context, 1)
        np.testing.assert_allclose(energy, expected_energy, rtol=rtol)
        np.testing.assert_allclose(forces, expected_forces, rtol=1e-5, atol=1e-6)

    def assert_same_settings(self, copy_, original):
        self.assertIsInstance(copy_, ColvarsCollectiveVariable)
        self.assertEqual(copy_.getConfigString(), original.getConfigString())
        self.assertEqual(copy_.getColvarName(), original.getColvarName())
        self.assertEqual(
            copy_.usesPeriodicBoundaryConditions(),
            original.usesPeriodicBoundaryConditions(),
        )
        self.assertEqual(copy_.getForceGroup(), original.getForceGroup())
        self.assertEqual(copy_.getName(), original.getName())

    def test_distance_value_and_forces(self):
        # Colvars reports distances in nm, as OpenMM does
        cv = ColvarsCollectiveVariable(colvar_config("d", DISTANCE), "d")
        reference = openmm.CustomBondForce("r")
        reference.addBond(0, 1, [])
        context = make_context(make_system(2, cv, reference), make_positions(2))
        self.assert_matches(context)

    def test_periodic_dihedral_value_and_forces(self):
        # Colvars reports dihedrals in degrees, wrapped into (-180, 180]; this is the
        # same convention as OpenMM's theta in radians
        for seed in [1, 2, 3]:
            with self.subTest(seed=seed):
                cv = ColvarsCollectiveVariable(colvar_config("phi", DIHEDRAL), "phi")
                reference = openmm.CustomTorsionForce(f"{RADIANS_TO_DEGREES}*theta")
                reference.addTorsion(0, 1, 2, 3, [])
                positions = make_positions(4, seed)
                context = make_context(make_system(4, cv, reference), positions)
                self.assert_matches(context)

    def test_orientation_angle_matches_applied_rotation(self):
        reference = make_positions(6)
        rotation, angle = rotation_matrix(1)
        center = reference.mean(axis=0)
        positions = (reference - center) @ rotation.T + center
        cv = ColvarsCollectiveVariable(orientation_angle(reference), "omega")
        context = make_context(make_system(6, cv), positions)
        energy, _ = energy_and_forces(context, 0)
        self.assertAlmostEqual(energy, angle, delta=1e-6)

    def test_getters_and_validation(self):
        cv = ColvarsCollectiveVariable(colvar_config("d", DISTANCE), "d")
        self.assertEqual(cv.getConfigString(), colvar_config("d", DISTANCE))
        self.assertEqual(cv.getColvarName(), "d")
        self.assertFalse(cv.usesPeriodicBoundaryConditions())
        cv.setUsesPeriodicBoundaryConditions(True)
        self.assertTrue(cv.usesPeriodicBoundaryConditions())
        with self.assertRaisesRegex(openmm.OpenMMException, "configString"):
            ColvarsCollectiveVariable("", "d")
        with self.assertRaisesRegex(openmm.OpenMMException, "colvarName"):
            ColvarsCollectiveVariable(colvar_config("d", DISTANCE), "")

    def test_serialization_round_trip(self):
        cv = configured_cv()
        xml = openmm.XmlSerializer.serialize(cv)
        self.assertIn("ColvarsCollectiveVariable", xml)
        # OpenMM returns plugin forces as plain Force objects, which work in a System
        restored = openmm.XmlSerializer.deserialize(xml)
        self.assertEqual(openmm.XmlSerializer.serialize(restored), xml)
        self.assert_same_settings(ColvarsCollectiveVariable._deserialize(xml), cv)
        positions = make_positions(4, 2)
        restored.setForceGroup(0)
        context = make_context(make_system(4, restored), positions)
        energy, forces = energy_and_forces(context, 0)
        expected = make_context(make_system(4, configured_cv()), positions)
        expected_energy, expected_forces = energy_and_forces(expected, 0)
        np.testing.assert_allclose(energy, expected_energy)
        np.testing.assert_allclose(forces, expected_forces)

    def test_python_copies(self):
        duplicates = {
            "copy": copy.copy,
            "deepcopy": copy.deepcopy,
            "pickle": lambda cv: pickle.loads(pickle.dumps(cv)),
        }
        for name, duplicate in duplicates.items():
            with self.subTest(name):
                cv = configured_cv()
                self.assert_same_settings(duplicate(cv), cv)

    def test_inside_a_custom_cv_force(self):
        expression = (
            "0.5*k*dphi^2 + d; dphi = phi - phi0 - 2*pi*floor((phi - phi0)/(2*pi) + 0.5)"
        ).replace("pi", str(math.pi))
        for platform in PLATFORMS:
            with self.subTest(platform=platform):
                distance = openmm.CustomBondForce("r")
                distance.addBond(0, 3, [])
                torsion = openmm.CustomTorsionForce("theta")
                torsion.addTorsion(0, 1, 2, 3, [])
                phi = ColvarsCollectiveVariable(DIHEDRAL_IN_RADIANS, "phi")
                coupling = coupling_force(phi, distance, expression)
                expected = coupling_force(torsion, copy.deepcopy(distance), expression)
                system = make_system(4, coupling, expected)
                context = make_context(system, make_positions(4, 3), platform)
                self.assert_matches(context)
                np.testing.assert_allclose(
                    coupling.getCollectiveVariableValues(context),
                    expected.getCollectiveVariableValues(context),
                    rtol=1e-6,
                )
                # New positions without a time step must give new values
                context.setPositions(make_positions(4, 4))
                self.assert_matches(context)

    def test_several_collective_variables_in_one_force(self):
        # Every instance owns its own Colvars module, so two of them must not interfere
        def torsion_cv(name, first):
            groups = "".join(
                f"    group{k + 1} {{ atomNumbers {atoms(first + k)} }}\n"
                for k in range(4)
            )
            component = "dihedral {\n    name t\n" + groups + "  }"
            config = colvar_config(
                name, component, "  customFunction t*0.017453292519943295\n"
            )
            return ColvarsCollectiveVariable(config, name)

        def torsion_reference(first):
            force = openmm.CustomTorsionForce("theta")
            force.addTorsion(first, first + 1, first + 2, first + 3, [])
            return force

        expression = "phi^2 + 2*psi^2"
        colvars = openmm.CustomCVForce(expression)
        expected = openmm.CustomCVForce(expression)
        for name, first in [("phi", 0), ("psi", 1)]:
            colvars.addCollectiveVariable(name, torsion_cv(name, first))
            expected.addCollectiveVariable(name, torsion_reference(first))
        context = make_context(make_system(5, colvars, expected), make_positions(5, 4))
        self.assert_matches(context)

    def test_invalid_configurations_are_rejected(self):
        cases = {
            "unknown colvar": (
                colvar_config("d", DISTANCE),
                "missing",
                "does not define a colvar",
            ),
            "bias": (
                colvar_config("d", DISTANCE)
                + "harmonic {\n  colvars d\n  centers 0.1\n  forceConstant 1.0\n}\n",
                "d",
                "must not define biases",
            ),
            "extended Lagrangian": (
                colvar_config(
                    "d",
                    DISTANCE,
                    "  extendedLagrangian on\n  extendedFluctuation 0.01\n"
                    "  extendedTimeConstant 200\n",
                ),
                "d",
                "extended Lagrangian",
            ),
            "vector": (
                colvar_config("v", DISTANCE.replace("distance {", "distanceVec {")),
                "v",
                "not a scalar",
            ),
            "scripted forces off": (
                "scriptedColvarForces off\n" + colvar_config("d", DISTANCE),
                "d",
                "scriptedColvarForces",
            ),
        }
        for case, (config, name, message) in cases.items():
            with self.subTest(case):
                cv = ColvarsCollectiveVariable(config, name)
                with self.assertRaisesRegex(openmm.OpenMMException, message):
                    make_context(make_system(2, cv), make_positions(2))

    def test_writes_no_files(self):
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as directory:
            os.chdir(directory)
            try:
                cv = ColvarsCollectiveVariable(colvar_config("phi", DIHEDRAL), "phi")
                context = make_context(make_system(4, cv), make_positions(4))
                context.getIntegrator().step(300)
                del context
                self.assertEqual(os.listdir(directory), [])
            finally:
                os.chdir(cwd)

    def test_user_log_level_is_accepted(self):
        # logLevel defaults to 0, and Colvars rejects it if set twice
        for config in ["logLevel 1\n", "  LOGLEVEL 1\n"]:
            with self.subTest(config=config):
                cv = ColvarsCollectiveVariable(config + colvar_config("d", DISTANCE), "d")
                reference = openmm.CustomBondForce("r")
                reference.addBond(0, 1, [])
                context = make_context(make_system(2, cv, reference), make_positions(2))
                self.assert_matches(context)

    def test_colvars_force_is_unchanged(self):
        config = colvar_config("d", DISTANCE) + (
            "harmonic {\n  name restraint\n  colvars d\n  centers 0.3\n"
            "  forceConstant 1000.0\n}\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            prefix = os.path.join(directory, "colvars")
            force = ColvarsForce.fromConfigString(config, prefix)
            positions = [openmm.Vec3(0, 0, 0), openmm.Vec3(0.35, 0, 0)]
            context = make_context(make_system(2, force), positions)
            energy, _ = energy_and_forces(context, 0)
            del context
        # 0.5 * 1000 kJ/mol/nm^2 * (0.35 nm - 0.30 nm)^2 = 1.25 kJ/mol
        self.assertAlmostEqual(energy, 1.25, places=6)


if __name__ == "__main__":
    unittest.main()
