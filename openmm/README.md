# OpenMM interface for the `chinfo` Colvars branch

This directory builds `openmm_colvars.ColvarsForce`, an ordinary OpenMM `Force`
that evaluates Colvars from the positions and periodic box of the current
OpenMM `Context`. Its bias energy is included in the potential energy and its
atomic forces are added to the forces used by the integrator. It provides
Reference/CPU implementations and a native CUDA bridge to this branch's
existing Colvars GPU atom groups and CV kernels.

The interface was adapted from upstream Colvars PR 979; provenance and license
information are in [UPSTREAM.md](UPSTREAM.md). The original `chinfo` snapshot did
not contain an OpenMM interface. This integration reuses the Colvars core and
`src/cuda` kernels already present in this repository.

## 1. Build and use from Python

The commands below run from the **Colvars repository root**. Use a Python
virtual environment or an activated conda environment so the build and runtime
use the same OpenMM installation. Local validation used OpenMM 8.6.1.

### Reference/CPU build

You need a C++17 compiler, CMake 3.24 or newer, and Python development headers. For a fresh Python
environment:

```bash
python -m venv .venv-openmm
source .venv-openmm/bin/activate
python -m pip install 'openmm==8.6.1' cmake ninja swig pytest

cmake -S openmm -B build-openmm -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPython_EXECUTABLE="$(command -v python)" \
  -DOPENMM_COLVARS_ENABLE_CUDA=OFF
cmake --build build-openmm --parallel 2
cmake --install build-openmm

python openmm/examples/minimal.py --platform Reference
python openmm/examples/minimal.py --platform CPU
```

CMake detects the headers and libraries shipped with the active pip OpenMM
installation, or the active conda prefix. If needed, pass
`-DOPENMM_ROOT=/path/to/openmm/prefix`, where that prefix contains `include/`
and `lib/`. To run directly from the build without installing:

```bash
PYTHONPATH="$PWD/build-openmm/python" python openmm/examples/minimal.py --platform Reference
```

The convenience script `openmm/scripts/build_install.sh` runs configuration,
build, and installation. It accepts extra CMake arguments and uses
`OPENMM_COLVARS_ENABLE_CUDA`, `OPENMM_COLVARS_BUILD_DIR`, and `BUILD_JOBS`.

### CUDA build

CUDA requires all of the following on the machine that will run the simulation:

- An NVIDIA GPU and a usable NVIDIA driver.
- An OpenMM installation containing its CUDA platform plugin and matching
  CUDA internal headers (`openmm/cuda/CudaContext.h`).
- A CUDA 12 or newer development toolkit with `nvcc`; the OpenMM CUDA runtime wheel alone
  does not supply the compiler needed to build the Colvars `.cu` sources.

For OpenMM 8.6.1, pip offers `cuda12` and `cuda13` extras. Choose the extra
compatible with the installed toolkit and driver. For example, with CUDA 12:

```bash
python -m pip install 'openmm[cuda12]==8.6.1'

cmake -S openmm -B build-openmm-cuda -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPython_EXECUTABLE="$(command -v python)" \
  -DOPENMM_COLVARS_ENABLE_CUDA=ON \
  -DCUDAToolkit_ROOT=/usr/local/cuda \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \
  -DCMAKE_CUDA_ARCHITECTURES=native
cmake --build build-openmm-cuda --parallel 2
cmake --install build-openmm-cuda

python openmm/examples/minimal.py --platform CUDA --precision mixed --colvars-mode gpu
```

Use `'openmm[cuda13]==8.6.1'` when building against a CUDA 13 toolkit instead.
For a compile-only build without a visible GPU, replace `native` with the
explicit compute capability of the intended GPU. Build architecture support
depends on the selected CUDA toolkit. An installed CPU-only interface rejects
CUDA context creation with a rebuild instruction; it does not silently download
the complete system into a generic host callback.

For `customFunction` expressions, provide the Lepton sources from the matching
OpenMM version:

```bash
git clone --depth 1 --branch 8.6.1 https://github.com/openmm/openmm.git ../openmm-source
cmake -S openmm -B build-openmm -DOPENMM_SOURCE_DIR="$PWD/../openmm-source"
cmake --build build-openmm --parallel 2
cmake --install build-openmm
```

Use the same source option with `build-openmm-cuda` for a CUDA build.

## 2. Minimal runnable example

[examples/minimal.py](examples/minimal.py) reads
[examples/colvars.in](examples/colvars.in), constructs two particles, and runs
1,000 steps of a 1 fs Verlet integration. It requires no PDB, topology, solvent,
or downloaded force-field files. The relevant configuration is:

```text
units gromacs
colvarsTrajFrequency 0
colvarsRestartFrequency 0

colvar {
  name distance
  width 0.1
  distance {
    group1 { atomNumbers 1 }
    group2 { atomNumbers 2 }
  }
}
harmonic {
  name restraint
  colvars distance
  centers 0.5
  forceConstant 1.0
}
```

Its Python setup is:

```python
from openmm import System, Context, Platform, VerletIntegrator
from openmm_colvars import ColvarsForce

system = System()
system.addParticle(12.0)
system.addParticle(24.0)
force = ColvarsForce("openmm/examples/colvars.in", "example")
force.setUseGpu(True)
force.setUsesPeriodicBoundaryConditions(False)
system.addForce(force)
integrator = VerletIntegrator(0.001)  # ps
platform = Platform.getPlatformByName("CUDA")
context = Context(system, integrator, platform, {"Precision": "mixed"})
context.setPositions([[0, 0, 0], [0.65, 0, 0]])  # nm
context.setVelocities([[0, 0, 0], [0, 0, 0]])
state = context.getState(getEnergy=True, getForces=True)
print(force.getBackend(context))
print(force.getColvarValue(context, "distance"))
integrator.step(1000)
context.getState(getEnergy=True)  # evaluate the final coordinates
force.flushOutputs(context)      # save state and pending Colvars output
```

Alternatively, use `ColvarsForce.fromConfigString(config_text, output_prefix)`.
Configure the force before constructing the `Context`. To use the compact
CPU Colvars path while retaining CUDA for OpenMM's physical forces, call
`force.setUseGpu(False)` before creating the context.

### Verify that GPU work is actually requested

The minimal example prints the **created Context's** platform, CUDA precision
and device properties, and the selected Colvars backend. With the device path
enabled the backend contains `CUDA/ColvarsGPU`; with `setUseGpu(False)` it
contains `CUDA/ColvarsCPU`. It also lists CV components that still execute on
the CPU. The existence of a registered `CUDA` platform alone is not evidence
that a simulation used that platform or that a CV used a GPU kernel.

Before integration, this example must report a CV of **0.65 nm**, bias energy
**1.125 kJ/mol**, and forces **(+15, 0, 0)** and **(-15, 0, 0) kJ/mol/nm**.
The final coordinates demonstrate that the bias participates in integration.
Use the CUDA tests and Nsight commands below to verify the device calculations
and data movement on an actual GPU.

## 3. Units, atom mapping, and periodic boundaries

| Quantity | Colvars configuration / interface | OpenMM |
| --- | --- | --- |
| Distance | nm (`units gromacs`) | nm |
| Bias energy | kJ/mol | kJ/mol |
| Atomic force | kJ/mol/nm | kJ/mol/nm |
| Mass | Da | Da |
| Atomic charge for charge-dependent CVs | Elementary charge (e), explicitly supplied | Particle order from the System |
| Temperature | K | K |
| Integration timestep | Colvars receives fs | Integrator uses ps; adapter multiplies by 1,000 |
| Atom number | One-based `atomNumbers` | Zero-based particle index |

There is no extra length, energy, or force scaling at the OpenMM boundary.
`units real` is rejected so a NAMD configuration cannot silently apply an
Å/kcal/mol restraint to nm/kJ/mol coordinates. Convert such configurations
explicitly, including CV widths, centers, and all bias parameters. PDB reference
coordinates use the PDB convention of Å and are converted to nm by the loader.

`dipoleMagnitude` and `dipoleAngle`, including uses inside composite CVs,
require `force.setParticleCharges(charges)` before constructing the context.
Supply exactly one finite numeric charge in elementary-charge units for every
System particle, in its zero-based particle order. These charges are fixed for
that context. The interface does not infer a charge model from OpenMM forces
or follow later changes to force parameters or global charge offsets. Missing
charges cause an initialization error for these components; distance and other
charge-independent CVs do not require a charge array. A dipole magnitude has
units e·nm in this interface.

For example, two particles with charges +1e and -1e use:

```python
force.setParticleCharges([1.0, -1.0])
force.setTargetTemperature(300.0)  # kelvin; default is 300 K
```

For a static `NonbondedForce` without parameter offsets, its particle charge
values can populate the array explicitly. With custom, polarizable, or
parameter-dependent interactions, supply the fixed effective-charge model
intended for the CV. Changing such interactions later does not update that
model. Set the Colvars target temperature to the simulation temperature when
using temperature-dependent features; it is configured separately from an
OpenMM thermostat.

For a scalar harmonic distance CV, Colvars uses

```text
k_effective = forceConstant / width^2
U = 0.5 * k_effective * (distance - center)^2
```

For groups A and B, a positive extension gives force
`k_effective*(r-r0)*(COM_B-COM_A)/r` on group A, and its negative on group B.
Each group force is distributed by atom mass fraction. The adapter adds the
resulting forces with the same sign to OpenMM's accumulator.

`PbcMode::Auto` enables PBC if another force in the same `System` declares
periodic boundary conditions. An isolated biased system can explicitly call
`setUsesPeriodicBoundaryConditions(True)`. Merely setting default box vectors
does not select periodic behavior in Auto mode. Each evaluation reads the
current box from the context, including box changes within one OpenMM step.
Both orthorhombic and OpenMM reduced triclinic boxes are covered by tests.

The CUDA bridge handles OpenMM's reordered particle indices, mixed-precision
position correction, and periodic cell offsets. It reconstructs the selected
coordinates with the same wrapping convention as `Context.getState` before
Colvars applies the component's minimum-image convention. As in other Colvars
engines, molecular groups must have a physically meaningful coordinate image;
the interface does not reconstruct a group arbitrarily split across images.

## 4. Repeated evaluation and supported biases

OpenMM may request energy, forces, or both several times at one step, including
at different positions or box vectors. Every requested evaluation recomputes
the CV and its gradients from the current coordinates. Bias state advances
only at the first evaluation for each OpenMM step count. Later evaluations
apply the same bias state to the current CV. They do not deposit additional
metadynamics hills, advance a moving restraint, or duplicate trajectory output.
Force-group queries that exclude `ColvarsForce` do not update it.

`getColvarValue(context, name)`, `getBackend(context)`,
`getBiasUpdateCount(context)`, and `getBiasState(context)` are diagnostics.
Call `context.getState(getEnergy=True)` or request forces first to evaluate the
current configuration; the diagnostics themselves do not advance bias state.
`getBiasUpdateCount` counts state-update evaluations, including a fresh run's
first evaluation at step zero, so it need not equal the number of integrated steps.
Colvars' default `stepZeroData off` does not deposit a metadynamics hill at
step zero. The analytical metadynamics tests explicitly use `stepZeroData on`
so their first hill has a known 2 kJ/mol height at the initial coordinates.

This integration supports the restraint variants and metadynamics whose state
updates are separated from force/energy evaluation. Unsupported bias types,
total-force sampling such as ABF/eABF, and `extendedLagrangian` CVs fail with an
explicit error. Engine-provided replica messaging, scripted forces, and a
multi-GPU CUDA context are unsupported. OpenCL and HIP are not exposed through this
OpenMM interface. Reference/CPU force group 31 is rejected because the OpenMM
`CustomCPPForce` implementation uses a signed shift; use groups 0–30 there.
CUDA uses an unsigned mask and accepts group 31.

Colvars' existing file-based `multipleReplicas` protocol remains available.
Local tests cover the overlap between a replaced replica state and its old
hill buffer, acceptance of genuinely new hills, and initially empty states.
Distributed multi-process replica runs have not been validated here.
The existing restrictions still apply: expanding grids and `keepHills on`
cannot be combined with `multipleReplicas`.

Colvars `timeStepFactor` must be 1 for every CV and bias. Larger factors use
Colvars-specific impulse scaling and sleeping-CV bookkeeping that this OpenMM
callback cannot honor for arbitrary energy/force queries; they are rejected
at initialization. Configure multiple time stepping through the OpenMM
integrator's force groups instead, with Colvars `timeStepFactor` left at 1.

Rewinding the step counter of a live context is rejected. Restart an adaptive
run by creating a new context with the required OpenMM state and a matching
Colvars restart supplied via `setRestartInputPrefix(...)`. OpenMM checkpoints
alone do not serialize this plugin's adaptive history.

The first evaluation after loading a Colvars state (including a state supplied
by `defaultInputStateFile`) represents the already evaluated saved step. It
recomputes forces without adding another hill or moving a restraint again.
Later OpenMM step increments continue from that saved Colvars step, even if
the recreated OpenMM context starts with a different step count. The first
restored geometry must agree with the saved CV within Colvars' restart check;
subsequent same-step coordinate changes are evaluated normally.

`getBiasState()` preserves the live metadynamics potential, including hills
that have not yet been projected onto grids. Such state files can contain an
optional `unprojectedHills` count; use this branch's matching build to resume
them exactly. The reader also accepts older state files without this field.
Do not assume a build predating this field can reproduce a pending-hill state.
Replica snapshots also record `lastHillStep`, the last deposited hill covered
by the snapshot, so a still-visible old buffer cannot add those hills again.
Old snapshots without this marker use the existing state-step cutoff convention.

Moving-restraint snapshots preserve current centers/walls/force constants,
their schedule, accumulated work, and partial-stage TI sums. Older moving-wall
states without `firstStep` are rejected because their schedule is ambiguous.
Older staged states without `dAdLambda` can resume at a stage boundary; a
mid-stage restart requires that saved sum. These errors prevent a restart
from silently changing the potential or losing accumulated statistics.

Call `force.flushOutputs(context)` after the final energy/force evaluation and
before destroying the context. It writes the current `.colvars.state` file,
PMF where configured, and pending hill-trajectory rows, including when a run
ends between output intervals or uses `outputFreq 0`. It preserves the bias
potential and update count; repeated calls and continuing integration are
supported. File errors are reported to Python at this explicit call. Context
destruction does not perform this finalization, so it must not substitute for
`flushOutputs()` when requesting output files.

## 5. GPU path and remaining host work

The adapter keeps persistent buffers for the **M selected atoms**, gathers their
coordinates on the device, and accumulates their forces into OpenMM's native
fixed-point force buffer on the device. It does not download or re-upload the
full N-particle coordinate or force arrays each step. The atom-order mapping
is refreshed after OpenMM reorder events; selected cell offsets are refreshed
when OpenMM changes them. Integration uses OpenMM's CUDA context and stream,
with events linking the existing Colvars per-group streams.

| Path / calculation | Host/device data movement per evaluation |
| --- | --- |
| CUDA with Colvars GPU atom groups / supported GPU CVs | Selected coordinates are gathered device-to-device; group/CV results required by the host bias are copied to host; small bias force coefficients return to the device. Atomic bias forces are accumulated on-device. |
| A CPU CV component within a GPU-enabled configuration | Existing Colvars fallback may copy that component's atom-group coordinates/gradients between host and device. The backend diagnostic lists CPU components. This is not an all-device CV calculation. |
| Distance CV in the minimal example and benchmark | Uses that existing CPU CV fallback with GPU atom groups; it is not a native GPU distance kernel. Native RMSD and coordination-number kernels in this branch use the device CV path. |
| CUDA with `setUseGpu(False)` | Downloads `3*M*sizeof(double) = 24*M` bytes of selected positions; force requests upload `24*M` bytes of selected forces. Energy-only requests do not upload atomic forces. |
| Bias algorithms | Bias state and bias energy/force coefficients are currently evaluated on the CPU, including harmonic restraints and metadynamics. |

Existing GPU group reductions require host waits before a CPU bias consumes
their results. Thus `CUDA/ColvarsGPU` denotes the device-coordinate and force
path, not a claim of a fully asynchronous, wholly device-resident bias. The
CUDA adapter allocates its buffers at context construction; allocation and
transfer behavior inside individual CVs must still be checked by profiling.

## 6. Numerical verification

Run the local regression suite after building:

```bash
PYTHONPATH="$PWD/build-openmm/python" python -m pytest -q -s openmm/tests \
  --openmm-platform Reference,CPU

PYTHONPATH="$PWD/build-openmm-cuda/python" python -m pytest -q -s openmm/tests \
  --openmm-platform CUDA --openmm-precision mixed --require-cuda
```

Repeat the CUDA command with `--openmm-precision double` and `single` when
validating those modes. `--require-cuda` makes an unavailable CUDA device a
failure. Without it, a requested but unavailable device is explicitly skipped;
the suite probes bare OpenMM first so a Colvars failure is never mistaken for
a missing GPU.

Tests cover unequal masses, nonconsecutive atom IDs, spectators, harmonic width
scaling, force direction and magnitude, Newton's third law, actual integration,
orthorhombic and triclinic minimum images, energy-only/force-only requests,
changed coordinates and box vectors at the same step, tiny coordinate changes,
force groups, metadynamics hill counts and fresh forces, and input rejection.
They also verify addition to physical forces, independent force instances,
overlapping context lifetimes and constructor failure, moving-restraint work,
well-tempered hill heights, step-counter rewind rejection, and coordination
pair-list invalidation after coordinate and box changes.
RMSD cases compare against independent Kabsch fits and finite-difference forces,
including nonzero fitting derivatives. Dipole cases check required explicit
charges, sparse charge/index mapping, and analytical bias energies/forces.
Analytical tests print actual values, expected values, backend, precision and
maximum force error as JSON when run with `-s`.

| Precision | Distance/RMSD CV tolerance (nm) | Energy tolerance (kJ/mol) | Force tolerance (kJ/mol/nm) |
| --- | ---: | ---: | ---: |
| Reference / CPU / CUDA double | 1e-12 | 1e-9 | 1e-8 |
| CUDA mixed | 2e-8 | 2e-6 | 2e-5 |
| CUDA single | 2e-6 | 2e-4 | 2e-3 |

These are absolute tolerances; CPU/Reference test precision is double for the
Colvars values and force calculations. GPU tolerances are test acceptance
criteria, not measured GPU errors. See [VALIDATION.md](VALIDATION.md) for what was
actually executed on the implementation machine.

## 7. GPU throughput and profiling

`examples/benchmark.py` builds a deterministic periodic Lennard-Jones system
with 32,768 atoms and two sparse groups of 128 atoms. It uses a 1 fs Verlet
integrator, a distance CV with a harmonic bias, identical starting coordinates
and velocities, fixed precision, and the same output frequency in every mode.
It measures a pure-OpenMM baseline, the Colvars GPU path, and the compact
Colvars CPU path **on the same CUDA platform**. It warms up before timing,
synchronizes at the boundaries, performs three measurements, and records GPU
properties, throughput, individual elapsed times and configuration as JSON.

```bash
python openmm/examples/benchmark.py \
  --platform CUDA --precision mixed --device 0 \
  --particles 32768 --group-size 128 \
  --warmup 1000 --steps 10000 --repeats 3 \
  --report-interval 0 --json benchmark-gpu.json
```

`--report-interval 0` means no reporting inside timed blocks, with an energy
read at their boundaries. To evaluate a reporting workload, repeat all three
modes with the same nonzero report interval. `ns/day = steps/s * dt_ps * 86.4`.
This microbenchmark does not predict throughput for a solvated protein or a
different CV. It deliberately fails if CUDA is unavailable.

Profile the modes separately, keeping all other arguments unchanged:

```bash
nsys profile --trace=cuda,nvtx,osrt --sample=none -o colvars-gpu \
  python openmm/examples/benchmark.py --modes gpu \
  --particles 32768 --group-size 128 --warmup 1000 --steps 10000 --repeats 1 \
  --json benchmark-profile-gpu.json

nsys profile --trace=cuda,nvtx,osrt --sample=none -o colvars-cpu \
  python openmm/examples/benchmark.py --modes cpu \
  --particles 32768 --group-size 128 --warmup 1000 --steps 10000 --repeats 1 \
  --json benchmark-profile-cpu.json

nsys stats --report cuda_api_sum,cuda_gpu_kern_sum,cuda_gpu_mem_time_sum colvars-gpu.nsys-rep
nsys stats --report cuda_api_sum,cuda_gpu_kern_sum,cuda_gpu_mem_time_sum colvars-cpu.nsys-rep
```

Inspect the steady-state interval after warmup in the Nsight timeline. Check
for full-system coordinate/force copies, per-step device allocation, selected
group transfers, and time spent in synchronization APIs. Compare the repeated
kernel and copy sequence with the table above. Kernel names include
`gatherColvarsPositions` and `addColvarsForces`, together with the reused
Colvars kernels. Include a `--modes baseline` profile when attributing waits
to OpenMM versus the Colvars adapter. Do not compare profiled throughput to
unprofiled throughput.

For a short CPU-only harness smoke test (not a GPU performance result):

```bash
python openmm/examples/benchmark.py --platform CPU --modes baseline,cpu \
  --particles 216 --group-size 4 --warmup 10 --steps 20 --repeats 1 \
  --json benchmark-cpu-smoke.json
```

No GitHub Actions or other GitHub CI is used by these commands.
