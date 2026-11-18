"""Fitted RMSD acceptance tests with an independent NumPy SVD oracle.

Coordinates/CV are nm, energies kJ/mol, and forces kJ/mol/nm.  The five-point
finite difference uses a 1e-5 nm displacement; comparisons use the shared
precision-specific engine tolerances (including 1e-8 for double forces).
"""

import json

import numpy as np
import openmm as mm
from openmm import unit
import pytest

from openmm_colvars import ColvarsForce


SELECTED = np.array([0, 2, 3, 5, 7])
SPECTATORS = np.array([1, 4, 6, 8])
REFERENCE = np.array([[0.10, 0.20, -0.30], [0.80, -0.10, 0.05],
                      [-0.25, 0.65, 0.20], [0.20, -0.35, 0.90],
                      [0.70, 0.55, 0.60]])
FIT_OFFSETS = np.array([[0.04, -0.01, 0.02], [-0.03, 0.02, 0.05],
                        [0.02, 0.06, -0.02], [-0.01, -0.04, 0.03],
                        [0.05, 0.01, -0.04]])
DEFORMATION = np.array([[0.03, -0.04, 0.02], [-0.05, 0.03, 0.01],
                        [0.04, 0.02, -0.03], [-0.02, -0.01, 0.06],
                        [0.01, 0.05, -0.04]])
WIDTH, CENTER, FORCE_CONSTANT = 0.20, 0.025, 4.0


def kabsch_rmsd(positions, fitting_reference):
    """Fit by SVD, then measure against the independently specified target."""
    selected = positions[SELECTED]
    centered = selected - selected.mean(axis=0)
    reference_center = fitting_reference.mean(axis=0)
    target = fitting_reference - reference_center
    left, _, right = np.linalg.svd(centered.T @ target)
    handedness = np.diag([1.0, 1.0, np.linalg.det(left @ right)])
    fitted = centered @ (left @ handedness @ right) + reference_center
    return float(np.sqrt(np.mean(np.sum((fitted - REFERENCE)**2, axis=1))))


def harmonic_energy(positions, fitting_reference):
    value = kabsch_rmsd(positions, fitting_reference)
    return 0.5 * FORCE_CONSTANT / WIDTH**2 * (value - CENTER)**2


def finite_difference_forces(positions, fitting_reference):
    forces = np.zeros_like(positions)
    step = 1.0e-5
    for atom in SELECTED:
        for direction in range(3):
            energies = []
            for multiple in (-2, -1, 1, 2):
                shifted = positions.copy()
                shifted[atom, direction] += multiple * step
                energies.append(harmonic_energy(shifted, fitting_reference))
            forces[atom, direction] = -(energies[0] - 8*energies[1]
                                        + 8*energies[2] - energies[3]) / (12*step)
    return forces


def coordinates(scale, angle, translation):
    # A nondegenerate, asymmetric five-atom geometry; neither rotation nor
    # translation alone suffices to align the deliberately deformed structure.
    axis = np.array([1.0, 2.0, -1.0]) / np.sqrt(6.0)
    cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
    rotation = (np.eye(3) * np.cos(angle)
                + (1-np.cos(angle)) * np.outer(axis, axis)
                + np.sin(angle) * cross)
    positions = np.zeros((9, 3))
    positions[SELECTED] = (REFERENCE + scale*DEFORMATION) @ rotation + translation
    positions[SPECTATORS] = [[2, 3, 4], [-1, 2, -3], [4, -3, 2], [3, 1, -2]]
    return positions


def inline_positions(values):
    return " ".join(f"({x:.16g}, {y:.16g}, {z:.16g})" for x, y, z in values)


@pytest.mark.parametrize("explicit_fit", [False, True],
                         ids=["minimum-rmsd", "nonzero-fitting-derivatives"])
def test_rmsd_harmonic_kabsch_and_atomic_forces(engine, tmp_path, explicit_fit):
    fitting_reference = REFERENCE + FIT_OFFSETS if explicit_fit else REFERENCE
    # Standard minimum RMSD disables fitting derivatives because they cancel.
    # A distinct fitting reference makes those derivatives nonzero and exercises
    # the fitting-gradient graph in addition to the native RMSD GPU component.
    fitting_options = ""
    if explicit_fit:
        fitting_options = f"""
      centerToReference on
      rotateToReference on
      enableFitGradients on
      refPositions {inline_positions(fitting_reference)}
"""
    config = f"""
units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {{
  name r
  width {WIDTH}
  rmsd {{
    atoms {{
      atomNumbers {" ".join(str(atom+1) for atom in SELECTED)}
      {fitting_options}
    }}
    refPositions {inline_positions(REFERENCE)}
  }}
}}
harmonic {{
  colvars r
  centers {CENTER}
  forceConstant {FORCE_CONSTANT}
}}
"""
    system = mm.System()
    # Unequal masses also distinguish the unweighted RMSD fit from a COM fit.
    for mass in (12, 7, 16, 24, 9, 14, 11, 32, 17):
        system.addParticle(mass)
    force = ColvarsForce.fromConfigString(config, str(tmp_path / "rmsd"))
    force.setUseGpu(engine.use_gpu)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName(engine.platform), engine.properties)
    maximum_errors = np.zeros(3)
    cv_tol, energy_tol, force_tol = engine.tolerances
    for scale, angle, translation in ((1.0, 0.57, [0.7, -0.4, 0.2]),
                                       (-0.8, -0.91, [-0.2, 0.8, -0.3]),
                                       (1.6, 1.13, [0.3, 0.1, 0.9])):
        positions = coordinates(scale, angle, translation)
        expected_cv = kabsch_rmsd(positions, fitting_reference)
        expected_energy = harmonic_energy(positions, fitting_reference)
        expected_forces = finite_difference_forces(positions, fitting_reference)
        context.setPositions(positions)
        # Repeat without stepping, then edit all coordinates in the next loop.
        # This catches stale fitted coordinates, gradients, or graph buffers.
        for _ in range(2):
            state = context.getState(getEnergy=True, getForces=True)
            value = force.getColvarValue(context, "r")
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            np.testing.assert_allclose(value, expected_cv, rtol=0, atol=cv_tol)
            np.testing.assert_allclose(energy, expected_energy, rtol=0, atol=energy_tol)
            np.testing.assert_allclose(forces, expected_forces, rtol=0, atol=force_tol)
            np.testing.assert_array_equal(forces[SPECTATORS], 0.0)
            np.testing.assert_allclose(forces.sum(axis=0), 0.0, rtol=0, atol=force_tol)
            maximum_errors = np.maximum(maximum_errors,
                [abs(value-expected_cv), abs(energy-expected_energy),
                 np.max(np.abs(forces-expected_forces))])
        assert force.getBiasUpdateCount(context) == 1
    backend = force.getBackend(context)
    if engine.platform == "CUDA" and engine.use_gpu:
        assert "CUDA/ColvarsGPU" in backend
        assert "CPU CVCs=none" in backend
        assert "selected atoms=5" in backend
    print(json.dumps({"platform": engine.platform, "gpu": engine.use_gpu,
                      "precision": engine.precision, "explicit_fit": explicit_fit,
                      "backend": backend, "maximum_absolute_errors": maximum_errors.tolist(),
                      "absolute_tolerances": [cv_tol, energy_tol, force_tol],
                      "finite_difference_step_nm": 1.0e-5}))
    del context, integrator
