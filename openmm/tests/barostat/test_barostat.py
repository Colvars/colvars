import filecmp
import gc
from pathlib import Path
import re
import sys
import traceback

import numpy as np
from openmm import (
    Context, CustomIntegrator, LangevinMiddleIntegrator, MonteCarloBarostat, NonbondedForce,
    Platform, System, Vec3,
)
from openmm.app import CharmmParameterSet, CharmmPsfFile, CutoffPeriodic, HBonds, PDBFile
from openmm.unit import (
    angstrom, bar, kelvin, kilojoule_per_mole, nanometer, picosecond, picoseconds,
)
from openmm_colvars import ColvarsForce


def make_context(config, prefix, barostat_freq, barostat_first=True,
                 colvars_group=0, no_barostat=False, custom_integrator=False):
    path = Path('../Common')
    psf = CharmmPsfFile(str(path / 'da.psf'))
    pos_pdb = PDBFile(str(path / 'da.coor.pdb'))
    vel_pdb = PDBFile(str(path / 'da.vel.pdb'))
    params = CharmmParameterSet(str(path / 'top_all22_prot.inp'),
                               str(path / 'par_all22_prot.inp'))
    psf.setBox(3.0*nanometer, 3.0*nanometer, 3.0*nanometer)
    system = psf.createSystem(params, nonbondedMethod=CutoffPeriodic,
                              nonbondedCutoff=1.0*nanometer, constraints=HBonds)
    for force in system.getForces():
        force.setForceGroup(0)
    barostat = MonteCarloBarostat(1*bar, 300*kelvin, barostat_freq)
    barostat.setRandomNumberSeed(5678)
    colvars = ColvarsForce(config, prefix)
    if custom_integrator:
        # Persist the last sample count and state for the twenty-step fixture.
        config_text = Path(config).read_text().replace('outputFreq 10', 'outputFreq 1')
        colvars.setConfigString(config_text.replace('colvarsRestartFrequency 50',
                                                   'colvarsRestartFrequency 1'))
    colvars.setRandomSeed(4321)
    colvars.setTargetTemperature(300)
    colvars.setForceGroup(colvars_group)
    if barostat_first and not no_barostat:
        system.addForce(barostat)
    system.addForce(colvars)
    if not barostat_first and not no_barostat:
        system.addForce(barostat)
    if custom_integrator:
        integrator = CustomIntegrator(0.001*picoseconds)
        integrator.addUpdateContextState()
        integrator.addComputePerDof('v', 'v+dt*f/m')
    else:
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond,
                                             0.001*picoseconds)
    integrator.setRandomNumberSeed(1234)
    context = Context(system, integrator, Platform.getPlatformByName('Reference'))
    context.setPositions(pos_pdb.positions)
    context.setVelocities(vel_pdb.positions * 0.1 * (nanometer/angstrom/picosecond))
    return context


def trajectory_steps(prefix):
    with open(prefix + '.colvars.traj') as stream:
        return [int(line.split()[0]) for line in stream
                if line.strip() and not line.startswith('#')]


def total_samples(prefix):
    with open(prefix + '.count') as stream:
        return sum(float(line.split()[-1]) for line in stream
                   if line.strip() and not line.startswith('#'))


def volume(context):
    box = context.getState().getPeriodicBoxVectors(asNumpy=True)
    return np.linalg.det(box.value_in_unit(nanometer))


def test_once_per_step():
    context = make_context('eabf.in', 'mc1', barostat_freq=1)
    initial_volume = volume(context)
    context.getIntegrator().step(200)
    final_volume = volume(context)
    del context
    gc.collect()
    assert trajectory_steps('mc1') == list(range(200))
    assert final_volume != initial_volume

    context = make_context('eabf.in', 'mc0', barostat_freq=0)
    context.getIntegrator().step(200)
    del context
    gc.collect()
    assert trajectory_steps('mc0') == list(range(200))
    assert total_samples('mc1') == total_samples('mc0')


def test_rejected_move():
    context = make_context('eabf.in', 'runA', barostat_freq=0)
    context.getIntegrator().step(1)
    for _ in range(99):
        context.getIntegrator().step(1)
    del context
    gc.collect()

    context = make_context('eabf.in', 'runB', barostat_freq=0)
    context.getIntegrator().step(1)
    for _ in range(99):
        st = context.getState(getPositions=True)
        box = st.getPeriodicBoxVectors()
        pos = st.getPositions(asNumpy=True)
        context.setPeriodicBoxVectors(*[v*1.01 for v in box])
        context.setPositions(pos*1.01)
        context.getState(getEnergy=True, getForces=True)
        context.setPeriodicBoxVectors(*box)
        context.setPositions(pos)
        context.getIntegrator().step(1)
    del context
    gc.collect()
    compared = 0
    for suffix in ('.colvars.traj', '.count', '.grad', '.colvars.state'):
        if Path('runA' + suffix).exists():
            assert filecmp.cmp('runA' + suffix, 'runB' + suffix, shallow=False), suffix
            compared += 1
    assert compared >= 3


def test_trial_energy():
    context = make_context('restraint.in', 'de', barostat_freq=0, colvars_group=1)
    context.getIntegrator().step(1)
    st = context.getState(getPositions=True)
    box = st.getPeriodicBoxVectors()
    pos = st.getPositions(asNumpy=True)
    for scale in (1.0, 1.01, 0.99):
        context.setPeriodicBoxVectors(*[v*scale for v in box])
        pos_s = pos*scale
        context.setPositions(pos_s)
        energy = context.getState(getEnergy=True, groups={1}).getPotentialEnergy()
        energy = energy.value_in_unit(kilojoule_per_mole)
        coords = pos_s.value_in_unit(nanometer)
        distance = np.linalg.norm(coords[4] - coords[0])
        assert abs(energy - 0.5*1000*(distance-0.30)**2) < 1e-6
        context.setPeriodicBoxVectors(*box)
        context.setPositions(pos)
    del context
    gc.collect()
    assert len(trajectory_steps('de')) == 1


def test_pre_step_queries_do_not_consume_rng():
    for prefix, queries in (('initial_no_query', 0), ('initial_queries', 3)):
        context = make_context('eabf.in', prefix, barostat_freq=0)
        for _ in range(queries):
            context.getState(getEnergy=True, getForces=True)
        context.getIntegrator().step(100)
        del context
        gc.collect()
    for suffix in ('.colvars.traj', '.count', '.grad', '.colvars.state'):
        assert filecmp.cmp('initial_no_query' + suffix, 'initial_queries' + suffix,
                           shallow=False), suffix


def make_extended_context(prefix, barostat_freq=0, barostat_seed=1, pressure=1,
                          restart_prefix=None, integrator=None):
    # Two separate molecules ensure scaling changes the CV, unlike an
    # intramolecular distance under scaleMoleculesAsRigid.
    system = System()
    for _ in range(2):
        system.addParticle(12.0)
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3))
    physical = NonbondedForce()
    physical.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
    physical.setCutoffDistance(1.0)
    for _ in range(2):
        physical.addParticle(0, 0.1, 0)
    system.addForce(physical)
    barostat = MonteCarloBarostat(pressure*bar, 300*kelvin, barostat_freq)
    barostat.setRandomNumberSeed(barostat_seed)
    system.addForce(barostat)
    colvars = ColvarsForce('unused')
    colvars.setConfigString('''
colvarsTrajFrequency 1
colvarsRestartFrequency 1
colvar {
  name d
  width 0.01
  extendedLagrangian on
  extendedFluctuation 0.01
  extendedTimeConstant 100
  extendedLangevinDamping 1
  outputAppliedForce on
  outputEnergy on
  outputVelocity on
  distance {
    group1 { atomNumbers 1 }
    group2 { atomNumbers 2 }
  }
}
''')
    colvars.setOutputPrefix(prefix)
    colvars.setRestartOutputPrefix(prefix)
    if restart_prefix is not None:
        colvars.setRestartInputPrefix(restart_prefix)
    colvars.setRandomSeed(4321)
    colvars.setTargetTemperature(300)
    colvars.setForceGroup(1)
    system.addForce(colvars)
    if integrator is None:
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
    integrator.setRandomNumberSeed(1234)
    context = Context(system, integrator, Platform.getPlatformByName('Reference'))
    context.setPositions([Vec3(0, 0, 0), Vec3(0.5, 0, 0)])
    context.setVelocities([Vec3(0, 0, 0), Vec3(0, 0, 0)])
    return context


def test_initial_extended_snapshot():
    context = make_extended_context('initial_spring')
    state = context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True)
    box = state.getPeriodicBoxVectors()
    # Obtain k without approximating the gas constant used by Colvars: the
    # two trial force measurements independently establish k and its constancy.
    measured_k = []
    for scale in (1.0, 1.01, 0.99, 1.01):
        context.setPeriodicBoxVectors(*[v*scale for v in box])
        context.setPositions(positions*scale)
        result = context.getState(getEnergy=True, getForces=True, groups={1})
        energy = result.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        force = result.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer)[0, 0]
        displacement = 0.5*(scale-1)
        if displacement == 0:
            assert abs(energy) < 1e-12 and abs(force) < 1e-12
        else:
            assert abs(energy - 0.5*force*displacement) < 1e-9
            measured_k.append(force/displacement)
    assert max(measured_k) - min(measured_k) < 1e-7
    assert measured_k[0] > 20000
    context.setPeriodicBoxVectors(*box)
    context.setPositions(positions)
    context.getIntegrator().step(1)
    del context
    gc.collect()
    assert trajectory_steps('initial_spring') == [0]

    # The user-approved boundary: a manually scaled first query defines q0.
    context = make_extended_context('initial_scaled_query')
    context.setPositions(positions*1.01)
    first = context.getState(getEnergy=True, groups={1})
    assert abs(first.getPotentialEnergy().value_in_unit(kilojoule_per_mole)) < 1e-12
    context.setPositions(positions)
    second = context.getState(getEnergy=True, groups={1})
    assert second.getPotentialEnergy().value_in_unit(kilojoule_per_mole) > 0.1
    context.getIntegrator().step(1)
    del context
    gc.collect()
    rows = [line.split() for line in Path('initial_scaled_query.colvars.traj').read_text().splitlines()
            if line.strip() and not line.startswith('#')]
    assert len(rows) == 1
    # Extended CV trajectories write actual d first, then auxiliary r_d.
    assert abs(float(rows[0][2]) - 0.505) < 1e-12


def test_first_mc_accept_and_reject():
    outcomes = set()
    # Opposite high pressures accept and reject the same fixed-seed proposal,
    # exercising real OpenMM paths without a private accept/reject test API.
    # Reference can share RNG state between the barostat and integrator, so
    # varying only the barostat seed need not vary its first proposal.
    for label, pressure in (('positive', 1000000), ('negative', -1000000)):
        prefix = 'initial_mc_' + label
        context = make_extended_context(prefix, barostat_freq=1,
                                        barostat_seed=1, pressure=pressure)
        before = volume(context)
        context.getIntegrator().step(1)
        after = volume(context)
        outcomes.add('accepted' if after != before else 'rejected')
        del context
        gc.collect()
        rows = [line.split() for line in Path(prefix + '.colvars.traj').read_text().splitlines()
                if line.strip() and not line.startswith('#')]
        assert len(rows) == 1 and int(rows[0][0]) == 0
        # Reported q is the pre-integration seed, not the accepted spatial d.
        assert abs(float(rows[0][2]) - 0.5) < 1e-12
    assert outcomes == {'accepted', 'rejected'}


def test_restart_before_first_mc():
    context = make_extended_context('initial_restart_source')
    # Normal restart output deliberately excludes step 0.
    context.getIntegrator().step(2)
    positions = context.getState(getPositions=True).getPositions()
    del context
    gc.collect()
    state = Path('initial_restart_source.colvars.state').read_text()
    q = float(re.search(r'\bextended_x\s+(\S+)', state).group(1))
    velocity = float(re.search(r'\bextended_v\s+(\S+)', state).group(1))
    assert velocity != 0
    context = make_extended_context('initial_restart_mc', barostat_freq=1,
                                    restart_prefix='initial_restart_source')
    context.setPositions(positions)
    context.getIntegrator().step(1)
    del context
    gc.collect()
    lines = Path('initial_restart_mc.colvars.traj').read_text().splitlines()
    labels = next(line.split()[1:] for line in lines if line.startswith('#'))
    rows = [line.split() for line in lines if line.strip() and not line.startswith('#')]
    assert len(rows) == 1
    assert abs(float(rows[0][labels.index('r_d')]) - q) < 1e-12
    assert abs(float(rows[0][labels.index('vr_d')]) - velocity) < 1e-12


def test_custom_integrator_force_cache():
    steps = 20
    observed = {}
    for label, queries in (('plain', 0), ('query', 3)):
        prefix = 'custom_cache_' + label
        context = make_context('eabf.in', prefix, barostat_freq=0,
                               custom_integrator=True)
        initial_positions = context.getState(getPositions=True).getPositions(asNumpy=True)
        for _ in range(queries):
            context.getState(getEnergy=True, getForces=True)
        context.getIntegrator().step(steps)
        final_positions = context.getState(getPositions=True).getPositions(asNumpy=True)
        assert np.array_equal(initial_positions.value_in_unit(nanometer),
                              final_positions.value_in_unit(nanometer))
        del context
        gc.collect()
        observed[label] = trajectory_steps(prefix)
    assert observed == {'plain': list(range(steps)), 'query': list(range(steps))}, observed
    assert total_samples('custom_cache_plain') == steps - 1
    assert total_samples('custom_cache_query') == steps - 1
    for suffix in ('.colvars.traj', '.count', '.grad', '.colvars.state'):
        assert filecmp.cmp('custom_cache_plain' + suffix, 'custom_cache_query' + suffix,
                           shallow=False), suffix


def test_trial_cache_consistency():
    integrator = CustomIntegrator(0.001*picoseconds)
    integrator.addPerDofVariable('saved', 0)
    names = ('committed', 'plus', 'minus', 'twice', 'cached')
    for name in names:
        integrator.addGlobalVariable(name, 0)
        integrator.addGlobalVariable('force_' + name, 0)
    integrator.addUpdateContextState()
    integrator.addComputePerDof('saved', 'x')
    integrator.addComputeGlobal('committed', 'energy1')
    integrator.addComputeSum('force_committed', 'f1*step(0.25-saved)')
    # Displacements exceed the position-cache tolerance. The final exact-input
    # trial must use the same auxiliary state as these uncached trials.
    for name, displacement in (('plus', 1e-5), ('minus', -1e-5),
                               ('twice', 2e-5), ('cached', 0)):
        integrator.addComputePerDof('x', 'saved+(%.15g)*step(saved-0.25)' % displacement)
        integrator.addComputeGlobal(name, 'energy1')
        integrator.addComputeSum('force_' + name, 'f1*step(0.25-saved)')
    context = make_extended_context('trial_cache_consistency', integrator=integrator)
    integrator.step(1)
    energies = {name: integrator.getGlobalVariableByName(name) for name in names}
    forces = {name: integrator.getGlobalVariableByName('force_' + name) for name in names}
    # Interpolate a quadratic from +d, -d, +2d at zero: weights 1, 1/3, -1/3.
    # The spring's force is linear, so averaging +/-d gives its zero value.
    energy_at_zero = energies['plus'] + (energies['minus'] - energies['twice'])/3
    force_at_zero = 0.5*(forces['plus'] + forces['minus'])
    del context
    gc.collect()
    assert trajectory_steps('trial_cache_consistency') == [0]
    assert abs(energies['cached'] - energy_at_zero) < 1e-10, (energies, energy_at_zero)
    assert abs(forces['cached'] - force_at_zero) < 1e-8, (forces, force_at_zero)


def test_startup_checks():
    try:
        context = make_context('restraint.in', 'o', 25, barostat_first=False)
    except Exception as error:
        assert 'after every Monte Carlo barostat' in str(error), str(error)
    else:
        del context
        gc.collect()
        raise AssertionError('Force ordering was not rejected')

    try:
        context = make_context('abf_plain.in', 'u', 25)
    except Exception as error:
        assert 'cannot be evaluated on trial configurations' in str(error), str(error)
    else:
        del context
        gc.collect()
        raise AssertionError('Plain ABF with a barostat was not rejected')

    context = make_context('abf_plain.in', 'no_mc', 25, no_barostat=True)
    context.getIntegrator().step(5)
    del context
    gc.collect()


if __name__ == '__main__':
    failed = False
    for test in (test_once_per_step, test_rejected_move, test_trial_energy,
                 test_pre_step_queries_do_not_consume_rng,
                 test_initial_extended_snapshot, test_first_mc_accept_and_reject,
                 test_restart_before_first_mc, test_custom_integrator_force_cache,
                 test_trial_cache_consistency, test_startup_checks):
        try:
            test()
        except Exception:
            failed = True
            print('FAIL:', test.__name__, flush=True)
            traceback.print_exc()
        else:
            print('PASS:', test.__name__, flush=True)
    sys.exit(1 if failed else 0)
