# Local validation record — 2026-11-18

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
