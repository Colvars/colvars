#!/usr/bin/env python3
"""Run an actual OpenMM integration with a Colvars distance restraint."""

import argparse
from pathlib import Path

import openmm as mm
from openmm import unit
from openmm_colvars import ColvarsForce


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("CUDA", "CPU", "Reference"), default="CUDA")
    parser.add_argument("--precision", choices=("single", "mixed", "double"), default="mixed")
    parser.add_argument("--colvars-mode", choices=("gpu", "cpu"), default="gpu")
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("colvars-minimal"))
    args = parser.parse_args()
    if args.steps < 0:
        parser.error("--steps must be nonnegative")

    system = mm.System()
    for mass in (12.0, 24.0):
        system.addParticle(mass)
    args.output.mkdir(parents=True, exist_ok=True)
    force = ColvarsForce(str(Path(__file__).with_name("colvars.in")), str(args.output / "colvars"))
    force.setUseGpu(args.colvars_mode == "gpu")
    force.setUsesPeriodicBoundaryConditions(False)
    force.setForceGroup(1)
    system.addForce(force)

    integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
    platform = mm.Platform.getPlatformByName(args.platform)
    properties = {"Precision": args.precision} if args.platform == "CUDA" else {}
    if args.platform == "CPU":
        properties["Threads"] = "1"
    context = mm.Context(system, integrator, platform, properties)
    context.setPositions([[0, 0, 0], [0.65, 0, 0]] * unit.nanometer)
    context.setVelocities([[0, 0, 0], [0, 0, 0]])
    state = context.getState(getEnergy=True, getForces=True)
    print("OpenMM version:", mm.__version__)
    print("OpenMM platform:", context.getPlatform().getName())
    print("Colvars backend:", force.getBackend(context))
    for name in ("Precision", "DeviceIndex", "DeviceName"):
        if name in platform.getPropertyNames():
            print(f"CUDA {name}:", platform.getPropertyValue(context, name))
    print("Initial distance: %.12g nm (expected 0.65 nm)" % force.getColvarValue(context, "distance"))
    print("Initial bias energy:", state.getPotentialEnergy(), "(expected 1.125 kJ/mol)")
    print("Initial forces:", state.getForces(), "(expected +15, -15 kJ/mol/nm along x)")

    integrator.step(args.steps)
    state = context.getState(getEnergy=True, getPositions=True)
    print("Integrated steps:", context.getStepCount())
    print("Bias state updates:", force.getBiasUpdateCount(context))
    print("Final distance:", force.getColvarValue(context, "distance"), "nm")
    print("Final potential energy:", state.getPotentialEnergy())
    print("Final positions:", state.getPositions())
    force.flushOutputs(context)
    print("Final Colvars state:", args.output / "colvars.colvars.state")


if __name__ == "__main__":
    main()
