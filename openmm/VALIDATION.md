# Local validation

## Windows hardware run, 2026-11-19

Starting revision `8379ea411672171af23544e3186cac2507c413d2`. Local Windows x64
builds used Python 3.11.18, OpenMM 8.6.1, MSVC 14.44, SDK 10.0.22621,
SWIG 4.5.1, Ninja 1.13.2, CMake 4.5.0 and NVCC 12.9.86. The CUDA device was
an NVIDIA GeForce RTX 5070 Ti (sm_120), driver 610.88. The OpenMM C++ headers,
import libraries, Python module and CUDA plugin came from the same conda
installation. Lepton came from OpenMM tag 8.6.1.

| Status | Check | Observation |
| --- | --- | --- |
| PASS | Separate CPU-only and CUDA builds and installs | Both compiled, linked, and imported from clean temporary installation directories in new Python processes. |
| PASS | Reference and CPU plugin suite | 190 passed. Includes real accepted and rejected MC proposals, adaptive/barostat exclusion, transactional box rejection, Unicode paths and protection of existing state after a failed backup. |
| PASS | CUDA plugin, GPU groups and compact CPU fallback | 186 passed in each of single, mixed and double. Includes device RMSD versus Kabsch/finite-difference forces, coordination numbers, restart and multiple force instances. |
| SKIP | Pathologically thin box | Four CUDA entries per precision (two box shapes × two backends). A height of 1e-180 nm exercises error handling on the host; successful rejection and byte-for-byte adaptive history preservation were checked on Reference and CPU. |
| PASS | Two-particle Reference and CPU example | Initial distance 0.65 nm, bias 1.1250000000000002 kJ/mol and x forces +15 and −15 kJ/mol/nm. Final distance after 1,000 steps 0.36159136586640184 nm. |
| PASS | Two-particle CUDA mixed integration | CUDA/ColvarsGPU reported GPU groups with CPU distance component. Initial bias 1.125000000000005 kJ/mol, x forces +15 and −15. Final distance 0.3615913658655028 nm, 1,001 updates including step zero. Trajectory and state written. |
| PASS | Adapter buffers on device | NVRTC-compiled gather and scatter ran on the GPU in single, mixed and double. Reordered indices, triclinic cell offsets, mixed corrections and accumulation into prefilled padded force arrays matched the oracle. Details in [adapter results](validation/rtx5070ti-adapter.json). |
| PASS | Colvars core with Lepton, GCC 16.2.0 (MinGW) | Six CTest unit tests and 108 functional executions. The new nearest-image test checked 5,000 displacements against an independent lattice enumeration (maximum error 1.65e-14 nm, tolerance 1e-12). |
| PASS | Golden outputs | 103 cases, 206 files, each compared with Spiff at relative tolerance 1e-6. |
| SKIP | Five functional comparisons | One disabled comparison and four cases without comparison files. All five functional executables ran successfully. Names are in the previous record below. |
| NOT TESTED | Distributed replica communication, multi-device simulations and full NAMD/GROMACS engine runs | File-based replica tests and shared-core tests were run locally. |
| NOT TESTED | Compute Sanitizer memcheck | The Windows driver debugger interface was disabled. The tool exited 77. The ordinary CUDA tests and two-context thread check ran. |

The first Windows replica run had 12 failures on missing temporary files.
After fixing file removal and replacement, all 12 passed and the full host
suite was rerun. The first CUDA link failed on the RMSD helper's mismatched
forward declaration. Both statuses above refer to the repaired builds.

The final plugin suites recorded 748 PASS, 12 SKIP and 0 FAIL. For the five
mass-weighted harmonic geometries (nonperiodic, orthogonal, general triclinic,
elongated tilted and tiny tilt), maximum absolute errors over both CUDA
backends were:

| Platform/precision | Distance (nm) | Bias energy (kJ/mol) | Atomic force (kJ/mol/nm) |
| --- | ---: | ---: | ---: |
| Reference or CPU | 0 | 0 | 1.42e-14 |
| CUDA single | 5.44e-8 | 2.12e-6 | 1.26e-4 |
| CUDA mixed | 4.44e-16 | 2.49e-14 | 2.33e-10 |
| CUDA double | 0 | 0 | 2.33e-10 |

The [oracle measurements](validation/hardware-oracles.json) include the
two-particle integration, PBC cases, MC acceptance/rejection and RMSD errors.
The [per-test results](validation/pytest-results.json) include the outcome and
platform for each case. For NumPy numerical assertions the table retains the
worst absolute discrepancy per source line and tolerance, the corresponding
actual and expected scalar, and the number of calls. Exception, state-text
and lifecycle assertions carry their test outcome.

Both plugin configurations used `cmake -S openmm -G Ninja`, Release,
`-DPython_EXECUTABLE=<env>/python.exe`, `-DOPENMM_ROOT=<env>/Library`, and
`-DOPENMM_SOURCE_DIR=<OpenMM 8.6.1 checkout>`. The CUDA build also used
`-DOPENMM_COLVARS_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120`. Commands
after setting `PYTHONPATH` to the respective build's `python` directory:

```powershell
python -m pytest -q -s openmm/tests --openmm-platform Reference,CPU
foreach ($precision in 'single','mixed','double') {
    python -m pytest -q -s openmm/tests --openmm-platform CUDA --openmm-precision $precision --require-cuda --results-json "results-$precision.json"
}
python openmm/examples/minimal.py --platform CUDA --precision mixed --colvars-mode gpu --steps 1000
```

Core configuration used `cmake -S cmake -G Ninja` with
`-DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ -DCMAKE_BUILD_TYPE=Release
-DBUILD_TOOLS=OFF -DBUILD_TESTS=ON -DBUILD_UNITTESTS=ON -DCOLVARS_OPENMP=OFF
-DCOLVARS_LEPTON=ON -DLEPTON_DIR=<OpenMM source>/libraries/lepton`.
After `cmake --build`, `ctest --output-on-failure` ran the six Windows
unit tests. The repository registers functional tests with CTest on Linux.
On Windows we copied `tests/input_files` to a temporary directory and ran
the same executable from that directory once for every non-Torch case:

```text
run_colvars_test --configuration_file CASE/test.in --trajectory_file trajectory.xyz --output_prefix CASE/test_out --force
spiff -r 1e-6 CASE/AutoDiff/FILE CASE/FILE
```

### RTX 5070 Ti throughput and transfer profile

One device, mixed precision, 32,768 particles, a 1 fs Verlet step and two
128-atom groups. Each mode used the same Lennard-Jones system, starting
positions, velocities, precision and reporting interval. The biased modes
added the same distance restraint. After 1,000 warmup steps, each timing
run measured three successive blocks of 10,000 steps. Overhead uses the
ratio of mean elapsed times. Profiling was a separate run.

| Mode | Mean steps/s before | Mean steps/s after | After range | After ns/day | After overhead |
| --- | ---: | ---: | ---: | ---: | ---: |
| OpenMM CUDA | 4,886.13 | 4,883.48 | 4,867.35–4,899.53 | 421.93 | 0% |
| CUDA/ColvarsGPU | 2,222.58 | 2,415.13 | 2,413.70–2,416.39 | 208.67 | +102.20% |
| CUDA/ColvarsCPU | 3,267.19 | 3,305.82 | 3,290.31–3,322.45 | 285.62 | +47.72% |

The device group reduction already supplies its center to a CPU distance
component. Omitting the redundant coordinate download increased GPU-mode
throughput by 8.66% on this system. Host bias evaluation, per-atom gradient
uploads and synchronization are still present. The compact CPU path completes
this particular distance workload faster. No throughput comparison was made
for RMSD or coordination numbers. See the [before timings](validation/rtx5070ti-before.json)
and [after timings](validation/rtx5070ti-after.json).

A separate three-mode run on the final installed DLL (after file I/O and
rejected-box repairs) measured **4,985.92 / 2,414.63 / 3,313.56 steps/s**
for baseline / GPU / compact CPU, or **430.78 / 208.62 / 286.29 ns/day**.
Overheads in that run were 106.49% and 50.47%. The records identify the
extension by SHA-256; the source checkout SHA identifies the parent of the
working changes. [Final timing samples](validation/rtx5070ti-release.json)

| Final trace, per integrated step | H2D bytes | D2H bytes | Executed kernels | Explicit host sync calls | Device alloc/free |
| --- | ---: | ---: | ---: | ---: | --- |
| OpenMM CUDA | 8,912.90 | 8,396.67 | 9.0010 | 1.0081 | 0 / 0 |
| CUDA/ColvarsGPU | 15,077.99 | 8,396.67 | 17.0018 | 3.0192 | 0 / 0 |
| CUDA/ColvarsCPU | 15,077.99 | 14,541.28 | 11.0012 | 2.0127 | 0 / 0 |

Nsight Systems 2026.5.2 captured the NVTX `integrate` range including the
final energy fence. Counts include OpenMM's periodic atom reorders. H2D and
D2H count explicit copy operations. GPU writes to mapped pinned host
memory (including group centers) do not appear as memcpy bytes. GPU-mode
D2H fell from 14,541.28 to 8,396.67 bytes/step. The difference is the selected
256 × 3 doubles plus the extra endpoint evaluation averaged across the run.
The JSON includes individual CUDA APIs, copy directions and kernel names.
[Before trace counts](validation/rtx5070ti-profile-before.json) ·
[After trace counts](validation/rtx5070ti-profile-after.json)

```powershell
python openmm/examples/benchmark.py --platform CUDA --precision mixed --device 0 --particles 32768 --group-size 128 --warmup 1000 --steps 10000 --repeats 3 --report-interval 0 --json timings.json
$env:NSYS_NVTX_PROFILER_REGISTER_ONLY = '0'
foreach ($mode in 'baseline','gpu','cpu') {
    nsys profile --trace=cuda,nvtx --cuda-graph-trace=node --sample=none --cpuctxsw=none --export=sqlite --capture-range=nvtx --capture-range-end=stop --nvtx-capture=integrate@openmm-colvars-benchmark -o "trace-$mode" python openmm/examples/benchmark.py --platform CUDA --precision mixed --device 0 --particles 32768 --group-size 128 --warmup 1000 --steps 10000 --repeats 1 --modes $mode --profile-range --json "profile-$mode.json"
}
```

## Previous local validation record, 2026-11-18

All builds and tests ran locally. GitHub Actions and other GitHub CI were not
used. Starting chinfo commit: 5b79a08438ad39413a0816ad5654fbd8f5ca57f8.

## Environment and verification boundaries

| Item | Observed environment |
| --- | --- |
| Python / OpenMM | 3.12.14 / 8.6.1 |
| C++ compiler / build | GNU 14.2.0, C++17, Release |
| Build tools | CMake 4.4.5, Ninja 1.13.2, SWIG 4.5.1 |
| Runtime platforms tested | Reference and CPU; CPU Threads=1 |
| Colvars arithmetic | cvm::real is double; OpenMM CPU has no selectable Precision property |
| NVIDIA GPU / driver | No /dev/nvidia devices, nvidia-smi, or libcuda.so.1 |
| Offline CUDA compilation | NVCC 13.0.88, runtime headers 13.0.96, architecture sm_80 |
| Embedded-kernel compilation | NVRTC 12.9.86, compute_80, single/mixed/double |
| Lepton / OpenMM source | Matching OpenMM tag 8.6.1 |

The CUDA Python extension compiled and linked against the installed OpenMM
8.6.1 CUDA library. This compile/link check used a CUDA 13 compiler for Colvars
and the available OpenMM CUDA 12 library. Runtime loading is blocked here by
the missing NVIDIA driver. For GPU execution, use matching OpenMM CUDA extras
and toolkit as documented in the README.

**No CUDA numerical result, GPU throughput, GPU model, or measured CUDA
transfer/synchronization profile is claimed here.** Compilation does not
replace the hardware acceptance tests below.

## Analytical distance and harmonic restraint

The six particles have masses [12, 11, 24, 36, 7, 48] Da. Groups use one-based
atom numbers 1 4 and 3 6: noncontiguous indexing and unequal masses. Particles
2 and 5 receive exactly zero Colvars force.

Configuration: units gromacs, distance width 0.2 nm, harmonic center 0.35 nm,
forceConstant 8.0 kJ/mol. Effective spring constant: 8 / 0.2² = 200 kJ/mol/nm².
The independent Python oracle computes mass-weighted centers, enumerates
nearby periodic images, evaluates the potential, and distributes equal and
opposite group forces by mass fraction.

Absolute tolerances: CV 1e-12 nm; energy 1e-9 kJ/mol; each atomic-force component
1e-8 kJ/mol/nm. Both Reference and CPU returned:

| Geometry | Measured CV (nm) | Measured energy (kJ/mol) | Maximum atomic-force error (kJ/mol/nm) |
| --- | ---: | ---: | ---: |
| Nonperiodic | 0.7348469228349533 | 14.810715401553251 | 1.4210854715202004e-14 |
| Orthorhombic boundary crossing | 0.22360679774997916 | 1.5975241575014667 | 3.552713678800501e-15 |
| Triclinic boundary crossing | 0.6403124237432851 | 8.428130337970073 | 1.7763568394002505e-15 |

CV and energy errors were zero at printed double precision.
[analytical_results.json](validation/analytical_results.json) contains the complete
input coordinates, masses, groups, box vectors, actual/expected CV, energy and
atomic-force arrays, errors and tolerances for all six cases. Tests also check
zero net bias force and individual force directions.

The two-particle example produced distance 0.65 nm, energy
1.1250000000000002 kJ/mol, and forces (+15,0,0) and (-15,0,0) kJ/mol/nm.
After 1,000 Verlet steps at 1 fs from zero velocities, the distance was
0.36159136586640184 nm. Separate one-step tests check coordinate changes
against the analytical Verlet update, including physical plus bias forces.
The final installed module was also used to run this CPU example for 1,000
steps; flushOutputs wrote a state file at step 1000 with that final CV.

## Test coverage

The OpenMM suite passed **164 tests** across Reference and CPU. Coverage:

- Analytical CV, energy, forces, PBC, file loading and actual integration.
- Physical-force addition, force-group exclusions, energy-only/force-only calls.
- Fresh results after same-step coordinate/box changes and small displacements.
- Ordinary and well-tempered metadynamics without duplicate hills.
- Analytical/gridded metadynamics restarts, legacy states, pending and boundary
  hills, read-only snapshots, and PMF output without live-state mutation.
- Restored-step queries and continuation through both the Python restart API
  and defaultInputStateFile, including independent OpenMM step numbering.
- Moving-restraint centers and accumulated work frozen during repeat evaluation.
- Moving walls and force-constant ramps across mid/end/staged restarts; preserved
  work, stage sums and nonzero schedule start; finite-difference wall TI derivatives.
- Two Colvars forces in one Context, independent overlapping Context lifetimes,
  and failed constructors while valid Contexts remain alive.
- Backward step-count rejection, invalid atom indices and unit mismatches.
- Coordination-number pairlists after arbitrary position/box changes.
- Dummy groups in either position, no selected atoms, and nonunit Colvars
  timeStepFactor rejection.
- File-based replica state/buffer replacement overlap, genuinely newer hills,
  initial empty snapshots and step-zero hills, without duplicate deposition.
- Explicit output flushing at default and nonzero output frequencies, with
  complete hill trajectories, restart/PMF files and unchanged adaptive state.
- Required charge provisioning for dipole CVs, correct explicit-charge forces,
  unequal masses, net charge, sparse reversed selections and zero-force spectators.
- RMSD with independent Kabsch fits and five-point finite-difference forces,
  including a distinct fitting reference and nonzero fitting derivatives.

Metadynamics tests explicitly set stepZeroData on; the library default stays off.
The core's **113 existing functional/unit CTest cases** also passed. These are
execution/smoke checks; the OpenMM analytical oracles supply explicit numerical
verification. In addition, the existing Spiff comparisons ran locally:
**103 cases compared 206 trajectory/force files**, all passing at the original
relative tolerance **1e-6**, with no relaxed-tolerance fallback. This includes
well-tempered metadynamics, grid expansion, rebinning and TI reference outputs.

CTest reports 108 successful Spiff entries, of which five intentionally skip
comparison under the existing scripts/reference layout: dipoleangle_harmonic-fixed,
distance_harmonicwalls-upper-moving, distance_harmonicwalls-upper-moving-staged,
eulerangles_ext_abf2d_euler, and eulerangles_ext_abf2d_euler_weighted. Spiff commit
e3ecd28540da9a495f36314a378f99c2e1e380cb was built locally with terminal highlighting
disabled because curses headers were unavailable. No golden files were changed.
No NAMD/GROMACS/LAMMPS engine simulation was run.
All core execution and golden comparisons were repeated successfully after
the final production changes, including the charge capability guard.

Independent review reproductions also confirmed the repairs directly: a replica
snapshot/buffer overlap now returns 6.0 kJ/mol instead of the erroneous 10.0;
a completed moving-wall ramp restores 4.5 kJ/mol and forces of magnitude
30 kJ/mol/nm, instead of reverting to the initial wall's 12.5 kJ/mol potential.

## Additional charge and fitting checks

For two 12 Da particles at 0 and 0.3 nm with explicit charges +1e/-1e,
dipoleMagnitude is 0.3 e·nm. A harmonic center of 0.2 e·nm, width of 0.1 e·nm
and forceConstant 1 gives 0.5 kJ/mol and forces +10/-10 kJ/mol/nm. Missing
charges now reject initialization, including when an OpenMM NonbondedForce
has charges. A sparse, unequal-mass, net-charged case also passed with CV
0.525 e·nm, energy 5.28125 kJ/mol and force magnitudes 56.875 kJ/mol/nm.

Both Reference and CPU RMSD tests use double precision, five sparse selected
atoms and four spectators, three transformed/deformed configurations, and
same-step repeat queries. The independent force oracle uses five-point
finite differences with h=1e-5 nm:

| RMSD case | Max CV error (nm) | Max energy error (kJ/mol) | Max force error (kJ/mol/nm) |
| --- | ---: | ---: | ---: |
| Standard minimum RMSD | 3.1919e-16 | 7.7022e-16 | 8.2133e-11 |
| Distinct fitting reference | 4.1633e-17 | 2.2204e-16 | 1.0695e-10 |
| Absolute tolerance | 1e-12 | 1e-9 | 1e-8 |

The fourth independent reviewer also ran two Reference contexts on two threads,
with 100 changing-coordinate/step evaluations per context. Values and update
counts remained independent; maximum harmonic energy error was 8.88e-16 kJ/mol.

## CUDA compile checks and inspected data path

The CUDA build compiled all five existing Colvars .cu translation units plus
the OpenMM adapter and linked the extension. Embedded gather/scatter source
also passed NVRTC compilation in single, mixed and double precision.

The adapter uses persistent selected-atom buffers, atom-reorder callbacks,
mixed-precision corrections, and selected cell offsets. It has no per-step
full-system coordinate/force download or device-buffer allocation. Existing
Colvars CV/atom-group code still transfers data and waits for host calculations,
as detailed in the README. Distance uses the existing CPU CV fallback with GPU
atom groups; existing GPU RMSD and coordination-number kernels are also wired
through this proxy.

A pinned-buffer race was corrected: host gradients are cleared only after the
already-required group event wait, so the preceding GPU force graph can finish
reading them. No extra barrier was added. This fix still needs hardware
execution and sanitizer validation.

Review also removed repeated construction of empty CUDA graphs for groups
without fitting gradients or active forces, and excluded dummy groups from
zero-sized atom-buffer graph operations. Constructor failures now destroy the
Colvars GPU module while the correct OpenMM CUDA context remains selected.
These paths compiled successfully; device execution is still outstanding.

A CPU smoke test of the benchmark used 216 Lennard-Jones particles, 10 warmup
steps and 20 timed steps. It verified the harness and JSON metadata only; these
short CPU timings are not GPU performance measurements.

## Reproduction and outstanding hardware tests

Build commands are in [README.md](README.md). Run from the repository root:

```bash
PYTHONPATH="$PWD/build-openmm/python" python -m pytest -q openmm/tests \
  --openmm-platform Reference,CPU

python openmm/scripts/check_cuda_sources.py \
  --openmm-source ../openmm-source --nvrtc-library /path/to/libnvrtc.so \
  --architecture compute_80

for precision in single mixed double; do
  PYTHONPATH="$PWD/build-openmm-cuda/python" python -m pytest -q openmm/tests \
    --openmm-platform CUDA --openmm-precision "$precision" --require-cuda
done

PYTHONPATH="$PWD/build-openmm-cuda/python" python openmm/examples/benchmark.py \
  --platform CUDA --precision mixed --modes baseline,gpu,cpu \
  --particles 32768 --group-size 128 --warmup 1000 --steps 10000 --repeats 3 \
  --report-interval 0 --json gpu-benchmark.json

PYTHONPATH="$PWD/build-openmm-cuda/python" nsys profile \
  --trace=cuda,nvtx,osrt --sample=none --force-overwrite=true -o colvars-gpu \
  python openmm/examples/benchmark.py --platform CUDA --precision mixed \
  --modes gpu --particles 32768 --group-size 128 \
  --warmup 1000 --steps 10000 --repeats 1 --json gpu-profile.json
nsys stats --report cuda_api_sum,cuda_gpu_kern_sum,cuda_gpu_mem_time_sum colvars-gpu.nsys-rep
```

Pending on hardware: single/mixed/double numerical agreement, GPU versus compact
CPU agreement, reorder/cell-offset correctness during long runs, stream/pinned
memory safety, actual GPU model and throughput for all three modes, and per-step
transfer/wait profiling. Use the identical system, integrator, precision, output
frequency, warmup and timed step counts in every mode. Repeat the profile with
--modes baseline and --modes cpu to isolate Colvars overhead.

To reproduce the existing core and golden-output checks, with the matching
OpenMM source already cloned as in the README:

```bash
cmake -S cmake -B ../build-colvars-core -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TOOLS=ON -DBUILD_TESTS=ON -DBUILD_UNITTESTS=ON \
  -DCOLVARS_LEPTON=ON -DLEPTON_DIR="$PWD/../openmm-source/libraries/lepton" \
  -DCOLVARS_OPENMP=OFF
cmake --build ../build-colvars-core --parallel 3
ctest --test-dir ../build-colvars-core --output-on-failure -E '_spiff$' --parallel 4

git clone --depth 1 https://github.com/Colvars/spiff.git ../spiff
make -C ../spiff TFLAG= TLIB=m
GIT_DIR="$PWD/.git" GIT_WORK_TREE="$PWD" PATH="$PWD/../spiff:$PATH" TERM=dumb \
  ctest --test-dir ../build-colvars-core --output-on-failure -R '_spiff$' --parallel 4
```
