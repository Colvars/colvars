# How to use the Colvars Dashboard in VMD

Note that some spiffy Dashboard features are only available in VMD compiled with
the Colvars Module version 2019-03-21 or later.

A more complete documentation can be found here:
http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:dashboard

## With VMD 1.9.4 or later (when available)

The plugin will be available under Extensions/Analysis.

## Adding the Dashboard to an older version of VMD

First, download the Colvars tree from https://github.com/Colvars/colvars/archive/master.zip

Note that `Makefile` is used within the VMD plugin distribution, by VMD's own build script.

Then there are two options:

### Installing on top of an existing VMD installation

Edit the PLUGINDIR variable in Makefile.local, and run
`make -f Makefile.local install`

To run it inside VMD, type in the VMD console:

```
package require cv_dashboard
cv_dashboard
```

If you don't have a personal `.vmdrc` file, copy the default `.vmdrc` file from
the directory where VMD is installed, into your home directory.
Then add the following lines at the end:

```
package require cv_dashboard
vmd_install_extension cv_dashboard cv_dashboard "Analysis/Colvars Dashboard"
```

The Dashboard will then be listed in the Extension/Analysis menu.

### Running directly from the Colvars source tree without installing

Type in the VMD console:
```
source <path/to/cv_dashboard.tcl>
cv_dashboard
```
