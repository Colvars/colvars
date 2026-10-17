from pathlib import Path

from openmm import (
    Context,
    HarmonicBondForce,
    Platform,
    System,
    Vec3,
    VerletIntegrator,
)
from openmm.unit import (
    dalton,
    kilojoule_per_mole,
    nanometer,
    picosecond,
)
from openmm_colvars import ColvarsForce

# Colvars writes trajectory/restart files using this prefix.
Path("output").mkdir(exist_ok=True)

system = System()
system.addParticle(39.9 * dalton)
system.addParticle(39.9 * dalton)

# Add physical forces before ColvarsForce.
bond = HarmonicBondForce()
bond.addBond(
    0,
    1,
    0.30 * nanometer,
    1000.0 * kilojoule_per_mole / nanometer**2,
)
system.addForce(bond)

# Colvars atomNumbers are one-based. Colvars length values use angstroms, so
# this 3.0 angstrom restraint corresponds to the 0.30 nm OpenMM bond length.
colvars_config = """
colvarsTrajFrequency 100

colvar {
  name d12
  width 0.1
  distance {
    group1 { atomNumbers 1 }
    group2 { atomNumbers 2 }
  }
}

harmonic {
  name restrain_d12
  colvars d12
  centers 3.0
  forceConstant 1.0
}
"""

colvars = ColvarsForce.fromConfigString(
    colvars_config,
    "output/colvars",
)
system.addForce(colvars)

integrator = VerletIntegrator(0.001 * picosecond)
platform = Platform.getPlatformByName("CPU")
context = Context(system, integrator, platform)
context.setPositions(
    [
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.3, 0.0, 0.0),
    ]
    * nanometer
)

integrator.step(1000)
state = context.getState(getEnergy=True, getPositions=True)
print("Platform:", context.getPlatform().getName())
print("Potential energy:", state.getPotentialEnergy())
print("Positions:", state.getPositions())
