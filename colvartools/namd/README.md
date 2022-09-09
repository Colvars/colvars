## Tools for replica-exchange computations in NAMD

The main script `colvars_replica_utils.tcl` can be sourced at any time during a NAMD job.  Inside it are defined Tcl functions that implement bias-exchange protocols; most of these functions rely on the Colvars scripting API.


### Examples

The folders below contain ready-to-run examples for each specific methodology.

Folder | Description
------ | -----------
`eq` | Standard MD run with multiple replicas; use to check NAMD build
`bemtd` | Example input for bias-exchange metadynamics
