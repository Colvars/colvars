"""Charge-dependent CVs require explicit charges in elementary-charge units (e).

The physical NonbondedForce is deliberately separate: it must not be mistaken
for a promise that Colvars has inferred a unique, fixed charge model from it.
Dipole magnitudes consequently have units e*nm in the native unit system.
"""

import numpy as np
import openmm as mm
import pytest

from openmm_colvars import ColvarsForce
from test_openmm_colvars import state_values


def configuration(component="dipoleMagnitude", selection="1 2"):
    if component == "dipoleMagnitude":
        block = f"dipoleMagnitude {{ atoms {{ atomNumbers {selection} }} }}"
    elif component == "dipoleAngle":
        block = """dipoleAngle {
    group1 { atomNumbers 1 2 }
    group2 { atomNumbers 1 }
    group3 { atomNumbers 2 }
  }"""
    else:
        raise ValueError(component)
    return f"""units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {{
  name dipole
  width 0.1
  {block}
}}
harmonic {{
  name restraint
  colvars dipole
  centers 0.2
  forceConstant 1
}}
"""


def make_context(engine, output, config, explicit_charges=None,
                 masses=(12.0, 12.0), physical_charges=(1.0, -1.0)):
    system = mm.System()
    physical = mm.NonbondedForce()
    physical.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    physical.setForceGroup(0)
    for mass, charge in zip(masses, physical_charges):
        system.addParticle(mass)
        physical.addParticle(charge, 0.3, 0.0)
    system.addForce(physical)
    force = ColvarsForce.fromConfigString(config, str(output))
    force.setUseGpu(engine.use_gpu)
    force.setUsesPeriodicBoundaryConditions(False)
    force.setForceGroup(1)
    if explicit_charges is not None:
        force.setParticleCharges(explicit_charges)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName(engine.platform), engine.properties)
    return context, integrator, force


@pytest.mark.parametrize("component", ("dipoleMagnitude", "dipoleAngle"))
def test_charge_dependent_cv_rejects_missing_explicit_charges(engine, tmp_path, component):
    # Even with a conventional +/-1e physical charge model, silently supplying
    # zero charges to Colvars would give a wrong CV and silently omit its force.
    with pytest.raises(Exception, match="setParticleCharges"):
        make_context(engine, tmp_path / "missing", configuration(component))


def test_explicit_charges_give_analytic_dipole_energy_and_forces(engine, tmp_path):
    context, integrator, force = make_context(
        engine, tmp_path / "dipole", configuration(), explicit_charges=[1.0, -1.0])
    context.setPositions([[0, 0, 0], [0.3, 0, 0]])
    # Only group1 is requested, isolating Colvars from the physical Coulomb
    # attraction in group0. The vector dipole is -0.3 e*nm along x; its magnitude
    # is 0.3 e*nm. k_eff=1/0.1^2 yields E=0.5 and bias forces +/-10.
    energy, forces = state_values(context, groups=1 << 1)
    np.testing.assert_allclose(force.getColvarValue(context, "dipole"), 0.3,
                               rtol=0, atol=engine.tolerances[0])
    np.testing.assert_allclose(energy, 0.5, rtol=0, atol=engine.tolerances[1])
    np.testing.assert_allclose(forces, [[10, 0, 0], [-10, 0, 0]],
                               rtol=0, atol=engine.tolerances[2])
    assert force.getBiasUpdateCount(context) == 1
    del context, integrator


@pytest.mark.parametrize("charges", ([1.0], [1.0, -1.0, 0.0]))
def test_explicit_charge_array_requires_all_system_particles(engine, tmp_path, charges):
    with pytest.raises(Exception, match="(?i)charge.*(?:one entry|particle|size|length)"):
        make_context(engine, tmp_path / "bad_length", configuration(), explicit_charges=charges)


def test_sparse_selection_maps_full_system_charges_and_mass_correction(engine, tmp_path):
    # Two selected particles occur at full-system indices0 and3; the requested
    # order is reversed and unrelated charged spectators sit between them.
    # The selected net charge is +1e, so d(p)/dx includes the mass-weighted COM
    # correction q_i-Q*m_i/M, not just the individual q_i.
    masses = (12.0, 5.0, 7.0, 36.0)
    charges = [2.0, 9.0, -7.0, -1.0]
    context, integrator, force = make_context(
        engine, tmp_path / "sparse", configuration(selection="4 1"),
        explicit_charges=charges, masses=masses, physical_charges=charges)
    context.setPositions([[0, 0, 0], [2, 1, 0], [-1, 2, 0], [0.3, 0, 0]])
    # COM_x=.225 nm and p_x=2*(-.225)-(.3-.225)=-.525 e*nm.
    # Effective charges are +/-1.75e, giving F_x=+/-56.875 kJ/mol/nm.
    energy, forces = state_values(context, groups=1 << 1)
    np.testing.assert_allclose(force.getColvarValue(context, "dipole"), 0.525,
                               rtol=0, atol=engine.tolerances[0])
    np.testing.assert_allclose(energy, 5.28125, rtol=0, atol=engine.tolerances[1])
    expected_forces = [[56.875, 0, 0], [0, 0, 0], [0, 0, 0], [-56.875, 0, 0]]
    np.testing.assert_allclose(forces, expected_forces, rtol=0, atol=engine.tolerances[2])
    np.testing.assert_array_equal(forces[[1, 2]], 0.0)
    del context, integrator
