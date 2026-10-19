from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, exit, stderr
from openmm_colvars import ColvarsForce


path = '../Common/'
psf = CharmmPsfFile(path + 'da.psf')
pos_pdb = PDBFile(path + 'da.coor.pdb')
vel_pdb = PDBFile(path + 'da.vel.pdb')
params = CharmmParameterSet(path + 'top_all22_prot.inp', path + 'par_all22_prot.inp')


system = psf.createSystem(params, nonbondedMethod=NoCutoff,
        nonbondedCutoff=1*nanometer, constraints=HBonds)

colvars = ColvarsForce("test.in", "test")
system.addForce(colvars)

#integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pos_pdb.positions)
simulation.context.setVelocities(vel_pdb.positions * 0.1 * (nanometer/angstrom/picosecond))

simulation.reporters.append(StateDataReporter(stdout, 1, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(20)
