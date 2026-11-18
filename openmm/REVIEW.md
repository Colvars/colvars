# Independent review and repair record

Review baseline: `5b79a08438ad39413a0816ad5654fbd8f5ca57f8` on `chinfo`.
Each numbered round reviews the complete current changes against that baseline,
including the OpenMM interface, core changes, build, examples, tests and docs.
Reviewers are separate agents; fixes and tests run locally. No GitHub CI is used.

## Round 1

The first independent complete review found issues requiring repair:

| Finding | Repair and local evidence |
| --- | --- |
| Analytical metadynamics restart discarded saved hills | Restore all explicit analytical hills; test both new and legacy state files. |
| Serializing gridded metadynamics changed live energy/forces | Preserve live grids and explicitly save the pending hill tail; verify snapshots and exact restart with keepHills off/on and boundary probes. |
| First restored evaluation could deposit a hill twice | Treat the saved step as already updated, then continue from it; test repeated and changed-coordinate queries and different OpenMM step numbering. |
| Nonunit Colvars timeStepFactor produced impulse-scaled forces and could abort teardown | Reject it on CVs and biases during initialization; test clean exception handling. |
| Failed CUDA construction could destroy Colvars GPU resources outside their context | Reset the proxy in a constructor-body catch while ContextSelector is alive; full CUDA compile succeeds. |

The review's suggestion concerning `defaultInputStateFile` was checked against
the parser: it was previously available only to GROMACS. This integration now
explicitly enables the existing option for OpenMM and tests it alongside
`setRestartInputPrefix()`.

Related repairs found while addressing the findings:

- Keep projected boundary hills and pending hills separate, avoiding duplicate
  off-grid energy/force contributions; preserve new list entries when replacing
  an existing history.
- Compute well-tempered hill heights from the complete current potential,
  including pending hills.
- Include pending hills in PMF output using a temporary grid, preserving live
  state and hill caches. Compare frequent output against a no-output control.
- Skip empty CUDA graphs for dummy groups, missing fitting gradients, and
  inactive force paths; avoid leaked graph handles and per-step graph creation.
- Recompute the inferred Lepton path on CMake reconfiguration, so adding
  OPENMM_SOURCE_DIR to an existing build enables customFunction support.

After these repairs, the complete OpenMM suite passed **74 tests** on Reference
and CPU. The existing core's 113 execution/unit cases passed; 103 golden-output
cases compared 206 trajectory/force files successfully at relative tolerance
1e-6 (five other comparison entries intentionally skip). All CUDA translation
units and the adapter extension compiled; GPU
execution could not be tested because this machine has no NVIDIA device/driver.
The detailed numerical results and hardware acceptance commands are in
[VALIDATION.md](VALIDATION.md).

## Round 2

A new independent reviewer inspected the complete updated patch and reproduced
three additional issues:

| Finding | Repair and local evidence |
| --- | --- |
| A replaced replica snapshot and its old buffer could count the same hills twice | Filter incremental records by an exact snapshot lastHillStep watermark; the original reproduction changed from 10.0 to the correct 6.0 kJ/mol. Tests also cover newer records and initially empty snapshots with legitimate step-zero hills. |
| Moving harmonicWalls restored their initial wall instead of the saved moving potential | Save and restore current walls, schedule, stage and accumulated work; a completed ramp restores 4.5 kJ/mol and force magnitude 30 kJ/mol/nm. Tests cover intermediate and completed ramps, both walls and force-constant ramps. |
| Final buffered hill output had no flush/finalization path | Add flushOutputs(context) through the core post_run output path; examples invoke it after the final evaluation. Tests verify complete hill rows and restart/PMF files without changing bias history or potential. |

Related fixes and safeguards:

- Access the last replica hill for logging only when parsing appended a hill.
- Initialize no-grid and replica-mirror fields explicitly; in particular,
  expand_grids was previously read uninitialized for no-grid replicas.
- Preserve the partial-stage dAdLambda sum. Reconstruct legacy moving walls
  only when their schedule is known, and reject legacy states lacking required
  schedule or mid-stage statistics instead of silently inventing them.
- Avoid indexing absent wall increments when accumulating work inside a wall.
- Correct the staged moving-wall TI derivative, checking its sign and width
  scaling against independent finite differences.

The coordinated local rebuild and complete OpenMM suite passed **148 tests**
on Reference/CPU. The updated CUDA library and SWIG module compiled and linked
successfully. The 113 core execution/unit tests and 103 actual golden-output
cases (206 files at relative tolerance 1e-6) were rerun successfully on these
sources, with the same five intentionally skipped comparison entries.
The independent CUDA source audit found no additional confirmed GPU defect;
the hardware limitations in VALIDATION.md remain.

## Round 3

A third independent reviewer inspected the full patch and reproduced one
additional correctness issue: dipole CVs silently used zero charges when the
caller omitted setParticleCharges, even when an OpenMM NonbondedForce had
nonzero charges. A two-particle case returned CV 0, energy 2 and zero force,
instead of CV 0.3 e·nm, energy 0.5 kJ/mol and forces of magnitude 10 kJ/mol/nm.

The repair adds an atomic-charge capability request to dipoleMagnitude and
dipoleAngle initialization, including composite uses. OpenMM requires an
explicit complete fixed-charge array; other engines retain their existing
charge behavior. Documentation states the units, particle order, fixed-charge
semantics and separate target-temperature setting. Twelve new Reference/CPU
cases verify the guard and analytical values, including sparse/net-charged
selections with unequal masses.

The review also identified a fitting-path coverage gap. Four Reference/CPU
RMSD cases now check independent Kabsch values and five-point finite-difference
forces for both standard RMSD and nonzero fitting derivatives. Maximum measured
force error was 1.0695e-10 kJ/mol/nm against a 1e-8 tolerance. No production
RMSD change was needed.

The final complete OpenMM suite passed **164 tests**. The core execution and
golden comparisons passed again, and the final CUDA extension compiled/linked.

## Round 4 — clean complete review

A fourth independent reviewer read the complete final production diff and
affected core context from the original baseline. A fresh adjunct independently
reviewed the build, bindings, examples, scripts, documentation, all tests, and
the numerical artifact. The final charge and RMSD tests were included.

**No new confirmed correctness, integration, compatibility or performance
blocker, and no additional required source fix, was found.** A separate
two-thread/two-Reference-context check ran 100 changed-coordinate/step evaluations
per context with independent CV values and update counts; maximum harmonic
energy error was 8.88e-16 kJ/mol. The complete diff passed whitespace checking.

This clean disposition is based on CPU execution and CUDA source/build review.
GPU runtime correctness, race detection, throughput and profiling still require
the hardware tests documented in VALIDATION.md. No GitHub CI was used in any
round; the final commit includes [skip ci] to suppress push-triggered workflows.
