## Bias-exchange metadynamics example

The `bemtd` folder contains an example input to run a bias-exchange metadynamics simulation using three collective variables.  These are the distances between the nitrogen atom of an ammonia molecule vs. the oxygen atom and the two hydrogen atoms of a water molecule.

The script `run.sh` will run two consecutive jobs, about an hour long each.  To change this, feel free to adjust the environment variable `numsteps` or range of values for the variable `ijob`.

### Analysis

Each replica runs continuously in Cartesian space, with the biasing potential switching from time to time.  This is *different* from what e.g. GROMACS+PLUMED uses, where the coordinates of the replicas are exchanged and represent trajectorioes with a constant bias instead.

To obtain the latter from the former, the `load_trajectories.py` script in the parent folder can be used to extract histograms and time series of selected variables with a constant bias, for example:
```
python3 ../load_trajectories.py --first 100000 --skip 10 --variables dist_N_O
```

Histogram files computed over trajectory segments without bias (e.g. `histogram.neutral-dist_N_O.dat`) can be compared directly with the equivalent files computed on simulations of the same systems without bias (example input in the `eq` folder).
