# Collective variables module (Colvars)

A software module for molecular simulation and analysis that provides a high-performance implementation of sampling algorithms defined on a reduced space of continuously differentiable functions (aka collective variables).

The module itself implements a variety of functions and algorithms, including free-energy estimators based on thermodynamic forces, non-equilibrium work and probability distributions.

## Obtaining and using

The easiest way to obtain binary versions of Colvars is via the simulation programs [NAMD](http://www.ks.uiuc.edu/Research/namd/) and [LAMMPS](http://lammps.sandia.gov/) and the visualization program [VMD](http://www.ks.uiuc.edu/Research/vmd/).  Please check [here](https://github.com/Colvars/colvars/releases) to see which version of Colvars is included with the round-number versions of VMD and NAMD.  Colvars is integrated with LAMMPS on a near-continuous basis, most often immediately after significant code changes.

## Documentation

The [Colvars webpage](http://colvars.github.io/) includes user documentation for the three codes, as well as a Doxygen-based [developer documentation](http://colvars.github.io/doxygen/html/).

The reference article is:
G. Fiorin, M. L. Klein, and J. Hénin, Molecular Physics 111, 3345 (2013).  
http://dx.doi.org/10.1080/00268976.2013.813594  \[[BibTex file](https://github.com/Colvars/colvars/blob/master/doc/ref_Fiorin_2013.bib?raw=true)\] \[[Endnote file](https://github.com/Colvars/colvars/blob/master/doc/ref_Fiorin_2013.ciw?raw=true)\]

## Example input

Colvars requires a configuration file, or alternatively configuration arguments given through scripting commands by the linked program.  In NAMD:
```
colvars on
cv configfile <Colvars configuration file>
```
In VMD (_Tip:_ try also the new "Colvars Dashboard" plugin):
```
cv molid top
cv configfile <Colvars configuration file>
```
In LAMMPS:
```
fix Colvars all colvars configfile <Colvars configuration file>
```
The contents of the configuration file are typically the same across all programs, for example:
```
colvar { # Define a new variable
  name d # Must give a name to this variable
  width 0.2 # Estimated fluctuation amplitude and/or grid resolution, "w_d"
  distance { # This variable is a distance between centers of mass (COMs)
    group1 { atomNumbers 1 2 3 } # List the atoms of the 1st group
    group2 { atomNumbers 4 5 6 } # List the atoms of the 2nd group
  }
}

harmonic { # Define a harmonic potential, 1/2*K*(d-d0)^2/w_d^2
  colvars d # Apply it to the variable "d"
  centers 5.0 # The center of the potential, "d0"
  forceConstant 10.0 # Force constant, "K"
}
```

Complete input decks for some of the most commonly used features are available in the `examples` repository:
https://github.com/Colvars/examples

See also the [examples](https://github.com/Colvars/colvars/tree/master/examples?raw=true) folder of this repository for other examples of configurations.  Configuration options (particularly, the selections of atoms) require minimal changes to reflect the specifics of each simulation.

The [tests](https://github.com/Colvars/colvars/tree/master/tests?raw=true) folder also contains functional segments of Colvars configuration, used to build numerical tests of code accuracy and stability.  Feel free to use these segments in your production runs.

## Updating to the latest version

To recompile each program with the most recent version of the module, [download](https://github.com/Colvars/colvars/archive/master.zip) the `master` branch of this repository, or clone it via git:
```
git clone https://github.com/Colvars/colvars.git
```
and run the provided `update-colvars-code.sh` script against the unpacked source tree of any of the supported programs:
```
./update-colvars-code.sh /path/to/NAMD_X.YY_Source ; # updates NAMD
./update-colvars-code.sh /path/to/vmd-X.Y.Z        ; # updates VMD
./update-colvars-code.sh /path/to/vmd-plugins      ; # updates VMD plugins
./update-colvars-code.sh /path/to/lammps           ; # updates LAMMPS
```
and recompile them.

The `update-colvars-code.sh` script and its supporting files are synchronized with the latest version of each program: [NAMD nightly build](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD), [VMD CVS](http://www.ks.uiuc.edu/Research/vmd/doxygen/cvsget.html) and the Github repository of [LAMMPS](https://github.com/lammps/lammps).  Earlier versions are not supported.

## Which version is recommended?

The `master` branch is to be considered the "*stable*" release at any given time; any bugfixes are released through `master` first.  The input syntax is near-completely *backward-compatible* and output files are *forward-compatible*.  Feel free to download Colvars and update NAMD, VMD or LAMMPS as needed.

Other branches are dedicated to the development of specific features: please use them at your own discretion.

## Feedback

Please use the "Issues" tab of this page to submit new bug reports or to suggest new features.

## License

This software is distributed under the GNU Lesser General Public License, version 3.  See COPYING.LESSER for complete licensing terms.
