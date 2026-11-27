"""Regression tests for the long-range correction of alchemical forces."""

import math
import os
import unittest

import openmm as mm
from openmm import unit
from openmm_fep import FEPAtoms, FEPConfig, prepare_fep_system, prepare_lambda_abf_system
from openmm_fep.openmm import set_context_lambda


PLATFORMS = os.environ.get("OPENMM_COLVARS_TEST_PLATFORMS", "Reference").split()
CONFIG = FEPConfig(
    alch_atoms=FEPAtoms(grow=[0, 1]),
    windows=[0.0, 0.25, 1.0],
    steps=[1, 1],
    alch_elec_lambda_start=0.0,
)


def two_particle_pme(dispersion_correction):
    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(3, 0, 0), mm.Vec3(0, 3, 0), mm.Vec3(0, 0, 3))
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    nonbonded.setPMEParameters(3.0, 24, 24, 24)
    nonbonded.setCutoffDistance(1.0)
    nonbonded.setUseDispersionCorrection(dispersion_correction)
    # Two nonexcluded charges in one alchemical partition generate the correction.
    for charge in (0.4, -0.4):
        system.addParticle(12.0)
        nonbonded.addParticle(charge, 0.28, 0.3)
    system.addForce(nonbonded)
    return system


class TestFEPLongRangeCorrection(unittest.TestCase):
    def check_system(self, prepare, lambdas, continuous=False):
        for dispersion in (True, False):
            prepared = prepare(two_particle_pme(dispersion), CONFIG)
            forces = prepared.system.getForces()
            custom = {
                force.getForceGroup(): force
                for force in forces
                if isinstance(force, mm.CustomNonbondedForce)
            }
            with self.subTest(dispersion=dispersion):
                # Check before Context creation so a regression reports a test failure
                # instead of aborting the process from OpenMM's LRC worker threads.
                correction = custom[prepared.force_group_electrostatics]
                self.assertFalse(correction.getUseLongRangeCorrection())
                self.assertEqual(custom[prepared.force_group_vdw].getUseLongRangeCorrection(), dispersion)
                [native] = [force for force in forces if isinstance(force, mm.NonbondedForce)]
                self.assertEqual(native.getUseDispersionCorrection(), dispersion)
                if continuous:
                    derivatives = {
                        correction.getEnergyParameterDerivativeName(i)
                        for i in range(correction.getNumEnergyParameterDerivatives())
                    }
                    self.assertIn(prepared.lambda_parameter_name, derivatives)
                for platform_name in PLATFORMS:
                    with self.subTest(platform=platform_name):
                        self.check_windows(prepared, platform_name, lambdas, continuous)

    def check_windows(self, prepared, platform_name, lambdas, continuous):
        platform = mm.Platform.getPlatformByName(platform_name)
        properties = {"Precision": "double"} if platform_name == "CUDA" else {}
        integrator = mm.VerletIntegrator(0.001)
        context = mm.Context(prepared.system, integrator, platform, properties)
        context.setPositions([mm.Vec3(0.3, 0.3, 0.3), mm.Vec3(0.9, 0.3, 0.3)])
        for lam in lambdas:
            with self.subTest(lam=lam):
                if continuous:
                    context.setParameter("lambda", lam)
                else:
                    set_context_lambda(context, prepared.config, lam)
                state = context.getState(energy=True, forces=True, parameterDerivatives=continuous)
                energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                forces = state.getForces().value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
                self.assertTrue(math.isfinite(energy))
                self.assertTrue(all(math.isfinite(x) for force in forces for x in force))
                if continuous:
                    self.assertTrue(math.isfinite(state.getEnergyParameterDerivatives()["lambda"]))

    def test_discrete_fep_midwindow(self):
        self.check_system(prepare_fep_system, (0.0, 0.25, 1.0))

    def test_continuous_lambda_midwindow_and_endpoint_derivatives(self):
        self.check_system(prepare_lambda_abf_system, (0.25, 1.0), continuous=True)


if __name__ == "__main__":
    unittest.main()
