#!/usr/bin/env python

# From a coordinate pdb and a velocity pdb, create a trr with a box dimension

import MDAnalysis

da = MDAnalysis.Universe("da.pdb")
vel = MDAnalysis.Universe("da.vel.pdb")

da.trajectory.ts.velocities = vel.trajectory.ts.positions
da.trajectory.ts.dimensions = [30, 30, 30, 90, 90, 90]

w = MDAnalysis.Writer("da.trr", da.trajectory.n_atoms)
w.write(da.trajectory)
w.close()
