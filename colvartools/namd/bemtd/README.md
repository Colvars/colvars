## Bias-exchange metadynamics example

The `bemtd` folder contains an example input to run a bias-exchange metadynamics simulation using three collective variables.  These are the distances between the nitrogen atom of an ammonia molecule vs. the oxygen atom and the two hydrogen atoms of a water molecule.

The script `run.sh` will run two consecutive jobs, about an hour long each.  To change this, feel free to adjust the environment variable `numsteps` or range of values for the variable `ijob`.

### Input requirements

The file `ammonia_water-bemtd.colvars.in` contains the Colvars configuration for all biases, including three metadynamics biases that can be exchanged, and one walls restraint that remains fixed.

Besides the required keywords, each metadynamics bias should contain also the following:
- `outputEnergy yes`, so that the biasing energy can be easily be read during post-processing of the `.colvars.traj` file under the column `E_<bias name>`;
- `outputFeatures active`, so that the active/inactive status of each bias is reported as 1/0 values in the `<bias name>_active` column; *when the value of this column is 0 for a given bias, its reported energy is out of date and should be ignored even if positive*.

The comments in the NAMD script `bemtd.namd` describe in more detail how to define the biases are subject to being exchanged.  The current implementation requires that *no more than one replica at a time is without biases*.


### Analysis

Each replica runs continuously in Cartesian space, with the biasing potential switching from time to time.  This is *different* from what e.g. GROMACS+PLUMED uses, where the coordinates of the replicas are exchanged and represent trajectorioes with a constant bias instead.

To obtain the latter from the former, the `load_trajectories.py` script in the parent folder can be used to extract histograms and time series of selected variables with a constant bias, for example:
```
python3 ../load_trajectories.py --first 100000 --skip 10 --variables dist_N_O --biases neutral
```

Histogram files computed over trajectory segments without bias (e.g. `histogram.neutral-dist_N_O.dat`) can be compared directly with the equivalent files computed on simulations of the same systems without bias (example input in the `eq` folder).
