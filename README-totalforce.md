### CHANGES IN THE DEFINITION OF SYSTEM FORCES (NOW TOTAL FORCES)

Starting from the version 2016-08-10 of the Colvars module, 
the role of system forces has been replaced by total forces.

These include *all* forces acting on a collective variable, whether they
come from the force field potential or from external terms
(e.g. restraints), including forces applied by Colvars itself.

In NAMD, forces applied by Colvars, IMD, SMD, TMD, symmetry
restraints and tclForces are now all counted in the total force.

In LAMMPS, forces applied by Colvars itself are now counted in the total
force (all forces from other fixes were being counted already).


### WHEN IS THIS CHANGE RELEVANT

This change affects results *only* when (1) outputSystemForce is
requested or (2) the ABF bias is used.  All other usage cases are
*unaffected* (colvar restraints, metadynamics, etc).

When system forces are reported (flag: outputSystemForce), their values
in the output may change, but the physical trajectory is never affected.
The physical results of ABF calculations may be affected in some cases.


### CHANGES TO ABF CALCULATIONS

Compared to previous Colvars versions, the ABF method will now attempt
to cancel external forces (for example, boundary walls) and it may be
not possible to resume through a state file a simulation that was
performed with a previous version.

There are three possible scenarios:

1. No external forces are applied to the atoms used by ABF: results are
unchanged.

2. Some of the atoms used by ABF experience external forces, but these
forces are not applied directly to the variables used by ABF
(e.g. another colvar that uses the same atoms, tclForces, etc): in this
case, we recommend beginning a new simulation.

3. External forces are applied to one or more of the colvars used by
ABF, but no other forces are applied to their atoms: you may use the
subtractAppliedForce keyword inside the corresponding colvars to
continue the previous simulation.
