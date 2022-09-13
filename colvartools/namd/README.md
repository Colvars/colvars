## Tools for replica-exchange computations in NAMD

The Tcl script `colvars_replica_utils.tcl` can be sourced at any time during a NAMD job: inside it are defined Tcl functions that implement bias-exchange protocols; most of these functions rely on the Colvars scripting API.


### Examples

The folders below contain ready-to-run examples for each specific methodology.

Folder | Description
------ | -----------
`eq` | Standard MD run with multiple replicas; use this example to test that the NAMD build being used supports replica-exchange protocols (e.g. multicore builds do not).
`bemtd` | Example input for bias-exchange metadynamics; requires an updated version of the Colvars module inside it (version >= xxxx-xx-xx).


### Post-processing


The Python script `load_trajectories.py` serves the purpose of reading `.colvars.traj` files produced during bias-exchange simulations (i.e. where the variables are continuous but the biases change), and to write out trajectories where the biases are constant.  This allows, for example to monitor the consistency of the statistics of each replica, in particular the neutral one (without biases).

Running this script in either the `eq` or `bemtd` folder should output the desired results.
