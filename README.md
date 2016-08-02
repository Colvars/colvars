Collective variables module (Colvars)
=======

A software module for molecular simulation and analysis programs that provides a flexible and high-performance platform for present and future algorithms.

The source code and binary are also distributed together with the simulation programs [NAMD](http://www.ks.uiuc.edu/Research/namd/) and [LAMMPS](http://lammps.sandia.gov/) and the visualization program [VMD](http://www.ks.uiuc.edu/Research/vmd/).

To recompile each program with the most recent version of the module, [download](https://github.com/colvars/colvars/archive/master.zip) the `master` branch of this repository, or clone it via git:
```
git clone https://github.com/colvars/colvars.git
```
and run the provided `update-colvars-code.sh` script against the unpacked source tree of each program:
```
./update-colvars-code.sh /path/to/NAMD_X.YY_Source ; # updates NAMD
./update-colvars-code.sh /path/to/vmd-X.Y.Z        ; # updates VMD
./update-colvars-code.sh /path/to/lammps-XXyyyZZ   ; # updates LAMMPS
```
The update script is synchronized with the latest version of each program: [NAMD nightly build](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD), [VMD CVS](http://www.ks.uiuc.edu/Research/vmd/doxygen/cvsget.html) and [LAMMPS development version](http://www.lammps.org/).

All bugfixes are released through the `master` branch, which is to be considered the "*stable*" release. 

Other branches are dedicated to the development of specific features: use them at your own discretion.

Please use the "Issues" tab of this page to submit new bug reports or to suggest new features.

This software is distributed under the GNU Lesser General Public License, version 3.  See COPYING.LESSER for complete licensing terms.

The webpage for this project is: http://colvars.github.io/

The reference article is:
G. Fiorin, M. L. Klein, and J. Hénin, Molecular Physics 111, 3345 (2013).  
http://dx.doi.org/10.1080/00268976.2013.813594  \[[BibTex file](https://github.com/colvars/colvars/blob/master/doc/ref_Fiorin_2013.bib?raw=true)\] \[[Endnote file](https://github.com/colvars/colvars/blob/master/doc/ref_Fiorin_2013.ciw?raw=true)\]
