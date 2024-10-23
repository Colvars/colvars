## Tools for replica-exchange computations in NAMD

The Tcl script `colvars_replica_utils.tcl` can be sourced at any time during a NAMD job: inside it are defined Tcl functions that implement bias-exchange protocols; most of these functions rely on the Colvars scripting API.

Unlike in the Tcl scripts in the `lib/replica` folder of a NAMD distribution tree, no `.history` files are written and un-shuffling the DCD trajectory files with `sortreplicas` is not yet supported.  This may change in the future.

In theory, any types of bias can be exchanged, but so far only metadynamics biases have been tested.

### Examples

The folders below contain ready-to-run examples for each specific methodology.  Running them requires a multiple-replica capable NAMD build; for a single workstation/node, this is usually a "netlrts" build (NB: NOT the "multicore" build).

Folder | Description
------ | -----------
`eq` | Standard MD run with multiple replicas; use this example to test that the NAMD build being used supports replica-exchange.
`bemtd` | Example input for bias-exchange metadynamics; requires an updated version of the Colvars module inside it (version >= xxxx-xx-xx).


### Post-processing

The Python script `load_trajectories.py` serves the purpose of reading `.colvars.traj` files produced during bias-exchange simulations (i.e. where the variables are continuous but the biases change), and to write out trajectories where the identity of the active bias in each replica is constant.  This allows, for example to monitor the statistical correctness of each replica with the given bias (or lack thereof).

Running this script in either the `eq` or `bemtd` folder should output the desired results.

Internally, this script uses the `plot_colvars_traj.py` script from the parent directory to read Colvars trajectory files.
