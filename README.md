Collective variables module
=======

A software module for molecular simulation programs, which provides a flexible and high-performance platform for present and future algorithms.

The source code is currently distributed together with the simulation programs NAMD and LAMMPS.

To obtain the latest version, download this repository (either as a Zip file or by using git), and run the update-colvars-code.sh script to update the program of your choice.

Use the "Issues" tab of this repository page to submit bug reports, or to suggest new features.



Fork Notes
==========
This fork has an implementation of linear restraints in additiona to harmonic.

Next up, create a subclass of the bias class which is similar to the
restraint class but which uses adapative linear. Remember, the force
constants update with SOGD and the gradient is approximated as the
running variance of the CV.