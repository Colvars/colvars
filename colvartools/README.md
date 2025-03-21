# Useful tools for Colvars simulations

This directory contains both standalone tools and Colvars scripts.

## Standalone tools
| File name | Summary |
| ------------- | ------------- |
| **abf_integrate** | Process free energy gradient from ABF/TI to generate a free energy surface using a legacy MCMC procedure. Superseded by Poisson integration (builtin or standalone) for dimensions 2 and 3, still needed for higher-dimension PMFs. Build using the provided **Makefile**.|
| **poisson_integrator** | Process free energy gradient from ABF/TI to generate a free energy surface as the solution of a Poisson equation. Build using the provided **Makefile**.|
| **poisson_integrator_conv** | Process free energy gradient from ABF/TI to generate a free energy surface as the solution of a Poisson equation, monitoring convergence towards a given discretized scalar field. Build using the provided **Makefile**.|
| **noe_to_colvars.py** | Parse an X-PLOR style list of assign commands for NOE restraints.|
| **plot_colvars_traj.py** | Select variables from a Colvars trajectory file and optionally plot them as a 1D graph as a function of time or of one of the variables.|
| **quaternion2rmatrix.tcl** | As the name says.|
| **test_scripted_gradients.tcl** | When implementing [colvars as scripted functions of components](http://colvars.github.io/colvars-refman-namd/colvars-refman-namd.html#sec:colvar_scripted), use this to test numerically the correctness of the analytical gradient. |
| **extract_weights_biases.py** | Script to read the weights and biases from a trained dense neural network model, and output them to plain text files suitable for the [NeuralNetwork CV](http://colvars.github.io/colvars-refman-namd/colvars-refman-namd.html#sec:neuralnetwork). Examples of how to train those models and save them on-the-fly can be found in the Supporting Information of [this paper](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.1c01010) |
## Colvars scripts

| File name | Summary |
| ------------- | ------------- |
| **abmd.tcl** | Adiabatic Biased MD after Marchi & Ballone JCP 1999 / Paci & Karplus JMB 1999; implemented in 20 lines of Tcl.|
| **pathCV.tcl** | Path-based collective variables after Branduardi et al. (2007). Optimized implementation that calculates only the necessary distances.|
