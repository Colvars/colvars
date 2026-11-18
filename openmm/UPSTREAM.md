# Upstream provenance

The initial OpenMM proxy, Force wrapper design, Python binding pattern, and PDB
selection code are adapted from the `openmm/` directory of Colvars pull request
979, commit 2ada654ecc7218409b66628f2a26c9deaba2f325.

The imported interface is distributed under GPL-3.0 (see LICENSE). The Colvars
core keeps its existing LGPL-3.0 license. This integration deliberately excludes
the unrelated FEP, TABF-integrator, and AlchemicalPME implementations.

## Source comparison, 2026-11-19

The GitHub API reported the same PR head on this date. The starting commit on
`Tinkbian/colvars:chinfo` was `8379ea411672171af23544e3186cac2507c413d2`.
The existing integration already contains the SWIG Force wrapper, file/string
configuration, explicit charges, state serialization, and native-unit PDB
loading. Both implementations convert PDB coordinates from Å to nm.

| Area | Official PR #979 | chinfo implementation and reason |
| --- | --- | --- |
| CUDA force callback | A deferred post-computation flushes the queue, reads host positions for the system, and uploads a full force array. | A direct ForceImpl callback gathers selected positions and scatters selected forces. The mapping follows OpenMM reorder events. |
| Repeated queries | A cache compares step, box, selected positions and alchemical parameters. | The adapter reevaluates the current geometry on each request. The core distinguishes state updates from evaluations on the saved bias surface. This supports changing a position or box at fixed step count. |
| Physical total forces | Deferred computation, ordering checks and group restrictions support sampling for certain configurations. | The selected-atom bridge runs before all physical forces are guaranteed to have completed. ABF requires ordering, group-completeness checks and unbiased total-force sampling on the correct step. The proxy rejects this request. |
| Extended dynamics and TABF | The PR supplies alchemical parameters and a dedicated TABF integrator. | Repeat-safe updates for the extra dynamical variables would require integrator support. Extended-Lagrangian CVs are rejected. |
| CV execution and bias | The PR's CUDA force path still evaluates Colvars using host coordinates. | GPU atom groups and native RMSD/coordination components reuse the branch's CUDA implementation. The distance component and host bias evaluation incur synchronization. See the backend diagnostic for a particular configuration. |

The PR description, implementation and branch head can be inspected at the
[upstream PR](https://github.com/Colvars/colvars/pull/979) and its
[pinned source](https://github.com/Colvars/colvars/tree/2ada654ecc7218409b66628f2a26c9deaba2f325/openmm).
