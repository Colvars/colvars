# Upstream provenance

The initial OpenMM proxy, Force wrapper design, Python binding pattern, and PDB
selection code are adapted from the `openmm/` directory of Colvars pull request
979, commit 2ada654ecc7218409b66628f2a26c9deaba2f325.

The imported interface is distributed under GPL-3.0 (see LICENSE). The Colvars
core keeps its existing LGPL-3.0 license. This integration deliberately excludes
the unrelated FEP, TABF-integrator, and AlchemicalPME implementations.
