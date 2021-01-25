# Collective variables module (Colvars)

A software module for molecular simulation and analysis that provides a high-performance implementation of sampling algorithms defined on a reduced space of continuously differentiable functions (aka collective variables).

First released in 2008 as part of the standard distribution of [NAMD](https://www.ks.uiuc.edu/Research/namd/) version 2.7b1, Colvars has also been integrated in [LAMMPS](https://lammps.sandia.gov/download.html) and [VMD](https://www.ks.uiuc.edu/Research/vmd/).

The module itself implements a variety of functions and algorithms, including free-energy estimators based on thermodynamic forces, non-equilibrium work and probability distributions.

## Obtaining and using

The easiest way to obtain pre-compiled versions of Colvars is via one of following:
- the molecular simulation program [LAMMPS](https://lammps.sandia.gov/download.html);
- the molecular simulation program [NAMD](https://www.ks.uiuc.edu/Research/namd/);
- the molecular visualization program [VMD](https://www.ks.uiuc.edu/Research/vmd/).

Please check [here](https://github.com/Colvars/colvars/releases) to see which version of Colvars is included with the round-number or "stable" versions of each code.

For the molecular simulation program [GROMACS](http://www.gromacs.org/), code may be compiled via our Colvars-patched [releases](#gromacs-colvars-releases).

## Documentation

The [Colvars webpage](https://colvars.github.io/) includes user documentation for the four codes, as well as a Doxygen-based [developer documentation](https://colvars.github.io/doxygen/html/).

The reference article is:
G. Fiorin, M. L. Klein, and J. Hénin, Molecular Physics 111, 3345 (2013).
https://dx.doi.org/10.1080/00268976.2013.813594  \[[BibTex file](https://github.com/Colvars/colvars/blob/master/doc/ref_Fiorin_2013.bib?raw=true)\] \[[Endnote file](https://github.com/Colvars/colvars/blob/master/doc/ref_Fiorin_2013.ciw?raw=true)\]

## Example input

Colvars requires a configuration file, or alternatively configuration arguments given through scripting commands by the linked program.
- In NAMD:
```
colvars on
cv configfile <Colvars configuration file>
```
- In VMD (_Tip:_ try also the new "Colvars Dashboard" plugin):
```
cv molid top
cv configfile <Colvars configuration file>
```
- In LAMMPS:
```
fix Colvars all colvars configfile <Colvars configuration file>
```
- In GROMACS:
```
gmx mdrun -s topol.tpr -deffnm topol -colvars <Colvars configuration file>
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
./update-colvars-code.sh /path/to/lammps           ; # updates LAMMPS
./update-colvars-code.sh /path/to/NAMD_X.YY_Source ; # updates NAMD
./update-colvars-code.sh /path/to/vmd-X.Y.Z        ; # updates VMD
./update-colvars-code.sh /path/to/vmd-plugins      ; # updates VMD plugins
./update-colvars-code.sh /path/to/gromacs-XXX.X    ; # update GROMACS
```
and recompile them.

The `update-colvars-code.sh` script support patching the latest version of each program:
- the [LAMMPS Github repository](https://github.com/lammps/lammps);
- the [NAMD Git repository](https://www.ks.uiuc.edu/Research/namd/development.html);
- the [CVS repositories of VMD and its plugins](https://www.ks.uiuc.edu/Research/vmd/doxygen/cvsget.html).

Earlier versions of these codes are not supported.

For GROMACS, pre-patched versions are provided below.

## Gromacs-Colvars releases

The following links allow to download several versions of Gromacs containing the latest version of Colvars:

 - **2020.5** in [Tar.gz](https://github.com/Colvars/gromacs/archive/v2020.5-colvars.tar.gz) and [Zip](https://github.com/Colvars/gromacs/archive/v2020.5-colvars.zip) formats

 - **2018.8** in [Tar.gz](https://github.com/Colvars/gromacs/archive/v2018.8-colvars.tar.gz) and [Zip](https://github.com/Colvars/gromacs/archive/v2018.8-colvars.zip) formats

Gromacs-Colvars versions not listed above are not supported.  When using the [Gromacs forum](https://gromacs.bioexcel.eu/) to discuss usage of these versions, please specify "GROMACS modification: **Yes**"

## Which version of Colvars is recommended?

The `master` branch is to be considered the "*stable*" release at any given time; any bugfixes are released through `master` first.  The input syntax is near-completely *backward-compatible* and output files are *forward-compatible*.  Feel free to download Colvars and update NAMD, VMD, LAMMPS or GROMACS as needed.

Other branches are dedicated to the development of specific features: please use them at your own discretion.

## Feedback

Please use the "Issues" tab of this page to submit new bug reports or to suggest new features.

## License

This software is distributed under the GNU Lesser General Public License, version 3.  See COPYING.LESSER for complete licensing terms.
