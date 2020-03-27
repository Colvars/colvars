# -*- tcl -*-

mgridforce yes

# Maps here are intentionally in a weird order (map2 before map, gridforces in
# the middle) to test the internal mappings

# Map computed from rotated coordinates

mgridforcepotfile       map2 ../Common/da-rotated.density.dx
mgridforcefile          map2 ../Common/da.atomicnumber.pdb
mgridforcecol           map2 B
mgridforcechargecol     map2 O
mgridforcescale         map2 0.0 0.0 0.0

# Fixed-scaling map used with GridForces (constant-scaling)

mgridforcepotfile       gridforces ../Common/gaussian.density.dx
mgridforcefile          gridforces ../Common/da.atomicnumber.pdb
mgridforcecol           gridforces B
mgridforcechargecol     gridforces O
mgridforcescale         gridforces 0.1 0.1 0.1

# Map computed from original coordinates

mgridforcepotfile       map ../Common/da.density.dx
mgridforcefile          map ../Common/da.atomicnumber.pdb
mgridforcecol           map B
mgridforcechargecol     map O
mgridforcescale         map 0.0 0.0 0.0
