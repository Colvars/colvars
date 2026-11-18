#!/usr/bin/env python3
"""Compare baseline, GPU Colvars and compact CPU Colvars on the SAME platform.

This is a reproducible periodic Lennard-Jones benchmark, not a prediction for
any protein system. CUDA is required by default; no silent CPU substitution.
"""

import argparse
import json
import math
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

import numpy as np
import openmm as mm
from openmm import unit
from openmm_colvars import ColvarsForce


def make_system(n, group_size):
    """Identical physical system and initial state for all benchmark modes."""
    if n < 2 * group_size:
        raise ValueError("--particles must be >= twice --group-size")
    side = math.ceil(n**(1.0 / 3.0))
    spacing = 0.45
    box = side * spacing
    if box <= 2.0:
        raise ValueError("At least 125 particles are required for a 1 nm periodic cutoff")
    grid = np.indices((side, side, side)).reshape(3, -1).T[:n] * spacing
    rng = np.random.default_rng(20261118)
    velocities = rng.normal(0.0, math.sqrt(0.008314462618 * 100.0 / 39.948), size=(n, 3))
    velocities -= velocities.mean(axis=0)
    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box, 0, 0), mm.Vec3(0, box, 0), mm.Vec3(0, 0, box))
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(1.0)
    nonbonded.setUseDispersionCorrection(False)
    for _ in range(n):
        system.addParticle(39.948)
        nonbonded.addParticle(0.0, 0.34, 0.997)
    system.addForce(nonbonded)
    # Sparse selections ensure the CUDA adapter cannot rely on identity ordering.
    first = np.arange(group_size) * 2
    second = np.arange(group_size) * 2 + 1
    delta = grid[second].mean(axis=0) - grid[first].mean(axis=0)
    delta -= box * np.rint(delta / box)
    center = float(np.linalg.norm(delta))
    config = f"""units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0
colvar {{
  name distance
  width 0.1
  distance {{
    group1 {{ atomNumbers {' '.join(str(i + 1) for i in first)} }}
    group2 {{ atomNumbers {' '.join(str(i + 1) for i in second)} }}
  }}
}}
harmonic {{
  name restraint
  colvars distance
  centers {center:.17g}
  forceConstant 0.01
}}
"""
    return system, grid, velocities, config


def synchronize(context, report):
    # The same explicit fence is used in every mode. Fetching energy guarantees
    # all previously submitted simulation work has completed before timing stops.
    state = context.getState(getEnergy=True)
    if report:
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(energy):
            raise RuntimeError("Benchmark produced non-finite potential energy")
    return state


def advance(integrator, context, steps, report_interval):
    if report_interval == 0:
        integrator.step(steps)
    else:
        remaining = steps
        while remaining:
            batch = min(remaining, report_interval)
            integrator.step(batch)
            synchronize(context, True)
            remaining -= batch


def device_metadata(platform, context):
    result = {"platform": platform.getName()}
    for name in platform.getPropertyNames():
        try:
            result[name] = platform.getPropertyValue(context, name)
        except mm.OpenMMException:
            pass
    if platform.getName() == "CUDA":
        try:
            result["nvidia_smi"] = subprocess.check_output([
                "nvidia-smi", "--query-gpu=index,name,uuid,driver_version,memory.total",
                "--format=csv,noheader"], text=True, timeout=15).strip().splitlines()
        except (OSError, subprocess.SubprocessError) as exc:
            result["nvidia_smi_unavailable"] = str(exc)
    return result


def run_mode(args, mode, output_dir):
    system, positions, velocities, config = make_system(args.particles, args.group_size)
    force = None
    if mode != "baseline":
        force = ColvarsForce.fromConfigString(config, str(output_dir / mode))
        force.setUseGpu(mode == "gpu")
        force.setUsesPeriodicBoundaryConditions(True)
        force.setForceGroup(1)
        system.addForce(force)
    integrator = mm.VerletIntegrator(args.timestep)
    platform = mm.Platform.getPlatformByName(args.platform)
    properties = {"Precision": args.precision, "DeterministicForces": "true"} if args.platform == "CUDA" else {}
    if args.platform == "CPU":
        properties["Threads"] = str(args.cpu_threads)
    if args.platform == "CUDA" and args.device is not None:
        properties["DeviceIndex"] = args.device
    context = mm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    context.setVelocities(velocities)
    advance(integrator, context, args.warmup, args.report_interval)
    synchronize(context, True)
    backend = force.getBackend(context) if force else "OpenMM only"
    if args.platform == "CUDA" and force:
        expected = "colvarsgpu" if mode == "gpu" else "colvarscpu"
        if expected not in backend.lower().replace(" ", ""):
            raise RuntimeError(f"Requested {mode} path, observed backend: {backend}")
    result = {"mode": mode, "backend": backend, "device": device_metadata(platform, context),
              "colvars_config": config if force else None,
              "samples_seconds": [], "samples_steps_per_second": []}
    for repeat in range(args.repeats):
        # Warmup has already happened. Endpoint fences are included equally for
        # all modes; no coordinate or force downloads occur unless Colvars needs
        # them internally. Report frequency is identical across modes.
        synchronize(context, True)
        start = time.perf_counter()
        advance(integrator, context, args.steps, args.report_interval)
        synchronize(context, True)
        elapsed = time.perf_counter() - start
        result["samples_seconds"].append(elapsed)
        result["samples_steps_per_second"].append(args.steps / elapsed)
    result["median_steps_per_second"] = statistics.median(result["samples_steps_per_second"])
    result["median_ns_per_day"] = result["median_steps_per_second"] * args.timestep * 86.4
    if force:
        result["bias_state_updates"] = force.getBiasUpdateCount(context)
        result["final_distance_nm"] = force.getColvarValue(context, "distance")
        # End-of-run output is outside every timed block in every mode.
        force.flushOutputs(context)
    del context, integrator
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("CUDA", "CPU", "Reference"), default="CUDA")
    parser.add_argument("--precision", choices=("single", "mixed", "double"), default="mixed")
    parser.add_argument("--modes", default="baseline,gpu,cpu", help="Comma-separated baseline,gpu,cpu")
    parser.add_argument("--particles", type=int, default=32768)
    parser.add_argument("--group-size", type=int, default=128)
    parser.add_argument("--timestep", type=float, default=0.001, help="Verlet timestep in ps")
    parser.add_argument("--warmup", type=int, default=1000)
    parser.add_argument("--steps", type=int, default=10000, help="Timed steps per repeat")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--report-interval", type=int, default=0, help="Energy read interval; 0 means endpoints only")
    parser.add_argument("--device", default=None, help="CUDA device index (single GPU)")
    parser.add_argument("--cpu-threads", type=int, default=1)
    parser.add_argument("--json", type=Path, default=Path("benchmark.json"))
    args = parser.parse_args()
    modes = args.modes.split(",")
    if not modes or any(mode not in ("baseline", "gpu", "cpu") for mode in modes):
        parser.error("--modes must contain baseline,gpu,cpu")
    if args.platform != "CUDA" and "gpu" in modes:
        parser.error("--modes gpu requires --platform CUDA; CPU smoke test: --modes baseline,cpu")
    if min(args.particles, args.group_size, args.steps, args.repeats, args.cpu_threads) <= 0 or args.timestep <= 0:
        parser.error("Particle/group sizes, steps/repeats/threads and timestep must be positive")
    if min(args.warmup, args.report_interval) < 0:
        parser.error("--warmup and --report-interval must be nonnegative")
    if args.device is not None and "," in args.device:
        parser.error("This benchmark tests a single CUDA device; --device accepts one index")

    report = {"openmm_version": mm.__version__, "numpy_version": np.__version__,
              "system": "periodic Lennard-Jones, 39.948 Da, sigma 0.34 nm, epsilon 0.997 kJ/mol, cutoff 1.0 nm",
              "integrator": "VerletIntegrator", "initial_velocity_seed": 20261118,
              "configuration": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              "results": []}
    with tempfile.TemporaryDirectory(prefix="colvars-openmm-benchmark-") as output:
        for mode in modes:
            result = run_mode(args, mode, Path(output))
            report["results"].append(result)
            print(json.dumps(result, indent=2), flush=True)
    baseline = next((r["median_steps_per_second"] for r in report["results"] if r["mode"] == "baseline"), None)
    if baseline:
        for result in report["results"]:
            result["throughput_relative_to_baseline"] = result["median_steps_per_second"] / baseline
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2) + "\n")
    print("Saved", args.json)


if __name__ == "__main__":
    main()
