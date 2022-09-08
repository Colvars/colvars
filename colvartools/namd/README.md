## Tools for replica-exchange computations in NAMD

The main script `colvars_replica_utils.tcl` can be sourced at any time during a NAMD job.  Inside it are defined Tcl functions that implement bias-exchange protocols; nearly all of these functions rely on the Colvars scripting API.


### Examples

The folders below contain ready-to-run examples for specific methodology.

Folder | Description
------ | -----------
`bemtd` | Example input for bias-exchange metadynamics
