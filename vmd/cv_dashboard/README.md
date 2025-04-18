When using the Dashboard, please cite the reference article:

J. Hénin, Laura J. S. Lopes, G. Fiorin, **Human learning for molecular simulations: the Collective Variables Dashboard in VMD**, J. Chem. Theo. Comput., 2022, 18, 3, 1945–1956 (https://pubs.acs.org/doi/10.1021/acs.jctc.1c01081).

# How to use the Colvars Dashboard in VMD

## Using the version included with VMD 1.9.4 alpha or later

Open the graphical interface by clicking the `Extensions/Analysis/Colvars Dashboard` menu item,
or by typing `cv_dashboard` in the VMD terminal.

Complete documentation is available here:
http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:dashboard

## Updating the Dashboard manually

The plugin can be updated manually to access recent features that are not yet included in official
releases of VMD.

**Using VMD 1.9.3 (released in 2016) is not recommended**.
as of 2021, VMD 1.9.4 "alpha" releases are the de facto current versions.

First, download the Colvars tree from https://github.com/Colvars/colvars/archive/master.zip
or clone the git repository.
Navigate to the `vmd/cv_dashboard` subdirectory.

### On a Unix-like system with GNU make

Edit the `DESTINATION` variable in Makefile.local with the path where VMD is installed,
and run `make -f Makefile.local install`
(`Makefile` is used within the VMD plugin distribution, by VMD's own build script, which is why we use
`Makefile.local`).
Then start VMD and open the interface using `Extensions/Analysis/Colvars Dashboard`.


### On Windows or if GNU make is not available
Find the path where VMD is installed
(by default: `C:\Program Files\University of Illinois\VMD`)
and replace the contents of `plugins/noarch/tcl/cv_dashboard1.x` with those of `vmd/cv_dashboard`
from the Colvars tree.
