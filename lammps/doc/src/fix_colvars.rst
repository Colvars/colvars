.. index:: fix colvars

fix colvars command
===================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID colvars *configfile* keyword value ...

* *ID*, *group-ID* are documented in :doc:`fix <fix>` command
* "colvars" = style name of this fix command
* *configfile* = configuration file for Colvars (use "*none*" to provide it inline)
* keyword = *output* or *input* or *unwrap* or *tstat* or *seed*

  .. parsed-literal::
     *output* value = state filename/prefix for Colvars (default: "out")
     *input* value = input state filename/prefix for Colvars (optional, default: "NULL")
     *unwrap* value = "yes" or "no" (default: "yes")
     *tstat* value = fix ID of thermostat applied to relevant atoms (default: "NULL")
     *seed* value = seed for random number generator (default: 1966)

Examples
""""""""

.. code-block:: LAMMPS

   # Create the fix using a config file, set prefix for its output files
   fix Colvars all colvars colvars.inp output ${JOB}

   # Communicate the LAMMPS target temperature to the Colvars module
   fix_modify Colvars tstat NPT

   # Add a new restraint specific to this LAMMPS run
   fix_modify Colvars config """
   harmonic {
       name restraint
       colvars distance1 distance2
       centers ${ref1} ${ref2}
       forceConstant ${kappa}
   }"""

Description
"""""""""""

This fix interfaces LAMMPS to the collective variables `Colvars
<https://colvars.github.io>`_ library, which allows to accelerate sampling of
rare events and the computation of free energy surfaces and potentials of
mean force (PMFs) for any set of collective variables using a variety of
sampling methods (e.g. umbrella-sampling, metadynamics, ABF...).

This documentation describes only the "fix colvars" command itself in
a LAMMPS script.  The Colvars library is fully documented in the included
`PDF manual <PDF/colvars-refman-lammps.pdf>`_ or in the webpage
`https://colvars.github.io/colvars-refman-lammps/colvars-refman-lammps.html
<https://colvars.github.io/colvars-refman-lammps/colvars-refman-lammps.html>`_.

The Colvars library is developed at `https://github.com/Colvars/colvars
<https://github.com/colvars/colvars>`_  A detailed discussion of its
implementation is in :ref:`(Fiorin) <Fiorin>`; additional citations for
specific features are printed at runtime if these features are used.

There are example scripts on the `Colvars website <https://colvars.github.io>`_
as well as in the ``examples/PACKAGES/colvars`` directory in the LAMMPS
source tree.

----------

The only required argument to the fix is the name of the Colvars
configuration file.  The contents of this file are independent from the MD
engine in which the Colvars library has been integrated, save for the units
that are specific to each engine.  In LAMMPS, the units used by Colvars are
consistent with those specificed by the :doc:`units <units>` command.

.. versionadded:: Colvars_2023-06-04 The special value "*none*"
                  (lowercase) initializes an empty Colvars module, which
                  allows loading configuration dynamically using
                  :doc:`fix_modify <fix_modify>` (see below).

The *group-ID* entry is ignored.  "fix colvars" will always apply to
the entire system, but specific atoms will be selected based on
selection keywords in the Colvars configuration file or files.  There is
no need to define multiple "fix colvars" instances and it is not
allowed.

The "output" keyword allows to specify the prefix of output files generated
by Colvars, for example "*output*.colvars.traj" or "output.pmf".  Supplying
an empty string suppresses any file output from Colvars to file, except for
data saved into the LAMMPS :doc:`binary restart <restart>` files.

The "input" keyword allows to specify an optional state file that contains
the restart information needed to continue a previous simulation state.
However, because "fix colvars" records its state in LAMMPS :doc:`binary
restart <restart>` files, this is usually not needed when using the
:doc:`read_restart <read_restart>` command.

The *unwrap* keyword controls whether wrapped or unwrapped coordinates are
passed to the Colvars library for calculation of the collective variables and
the resulting forces.  The default is *yes*, i.e. the image flags are used to
reconstruct the absolute atom positions.  Setting this to *no* will use the
current local coordinates that are wrapped back into the simulation cell at
each re-neighboring step instead.  For information about when and how this
affects results, please see
`https://colvars.github.io/colvars-refman-lammps/colvars-refman-lammps.html#sec:colvar_atom_groups_wrapping
<https://colvars.github.io/colvars-refman-lammps/colvars-refman-lammps.html#sec:colvar_atom_groups_wrapping>`_.

The *tstat* keyword can be either "NULL" or the label of a thermostatting
fix that thermostats all atoms in the fix colvars group. This will be
used to provide the colvars module with the current thermostat target
temperature.

The *seed* keyword contains the seed for the random number generator
that will be used in the colvars module.


Restarting
""""""""""

This fix writes the current state of the Colvars module into :doc:`binary
restart files <restart>`.  This is in addition to the text-mode
".colvars.state" state file that is written by the Colvars module itself.
The information contained in both files is identical, and the binary LAMMPS
restart file is also used by fix colvars when :doc:`read_restart
<read_restart>` is called in a LAMMPS script.  In that case, there is
typically no need to specify the *input* keyword.

As long as LAMMPS binary restarts will be used to continue a simulation, it
is safe to delete the ".colvars.state" files to save space.  However, when a
LAMMPS simulation is restarted using :doc:`read_data <read_data>`, the
Colvars state file must be available and loaded via the "input" keyword or
via a "fix_modify Colvars load" command (see below).

When restarting, the fix and the Colvars module should be created and
configured using either the original configuration file(s).


Output
""""""

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the energy due to all
external potentials defined in the Colvars configuration.  The scalar value
calculated by this fix is "extensive".

Aside from the state information in a ".colvars.state" file, other
`output files <https://colvars.github.io/colvars-refman-lammps/colvars-refman-lammps.html#sec:colvars_output>`_
are produced by Colvars depending on the type of simulation.
For this reason, the "output" keyword is required for fix colvars.


Controlling Colvars via `fix_modify`
""""""""""""""""""""""""""""""""""""

The :doc:`fix_modify <fix_modify>` command may be used on "fix colvars" in
either one of two ways:

(1) Provide updated values for the fix parameters, such as *output*, *input*,
    *unwrap*, *tstat* and *seed*.  Additionally, the :doc:`fix_modify
    <fix_modify>` *energy* keyword is supported by this fix to add the energy
    change from the biasing force added by Colvars to the global potential
    energy of the system as part of :doc:`thermodynamic output <thermo_style>`
    (the default is :doc:`fix_modify energy no <fix_modify>`).
    For example, in a multi-step LAMMPS script involving multiple thermostats
    (e.g. fix nvt followed by fix npt), Colvars can read a new thermostat's
    target temperature like this:

   .. code-block:: LAMMPS

      fix NVT all nvt ...
      fix Colvars all colvars <configfile> output equil1 tstat NVT
      run <NUMSTEPS>
      unfix nvt
      fix NPT all n ...
      fix_modify Colvars tstat NPT
      fix_modify Colvars output equil2


(2) .. versionadded:: Colvars_2023-06-04 Call one of the scripting
    functions provided by the Colvars module itself (a full list is available
    in the Colvars doc).  The arguments to these functions are provided as
    strings and passed to Colvars.

    LAMMPS variables referenced by their string representation
    "${variable}" will be expanded immediately.  Note also that this
    variable expansion *will also happen within quotes*, similar to what the
    :doc:`mdi <mdi>` command provides.  This feature makes it possible to use
    the values of certain LAMMPS variables in Colvars configuration strings.
    For example, to synchronize the LAMMPS and Colvars dump frequencies:

   .. code-block:: LAMMPS

      variable freq index 10000
      dump myDump all atom/zstd ${freq} dump.atom.zstd
      fix_modify Colvars config "colvarsTrajFrequency ${freq}"

.. note::

   Although it is possible to use :doc:`fix_modify <fix_modify>` at any time,
   its results will only reflect the state of the Colvars module at the end
   of the most recent "run" or "minimize" command.  Any new configuration
   added via "fix_modify Colvars configfile" or "fix_modify Colvars config"
   will only be loaded when the simulation resumes.  Configuration files or
   strings will be parsed in the same sequence as they were provided in the
   LAMMPS script.


Restrictions
""""""""""""

This fix is provided by the COLVARS package and is only available if LAMMPS
was built with that package (default in most builds).  Some of the features
also require code available from the LEPTON package.  See the :doc:`Build
package <Build_package>` page for more info.

There can only be one Colvars instance defined at a time.  Since the
interface communicates only the minimum required amount of information, and
the Colvars module itself can handle an arbitrary number of collective
variables, this is not a limitation of functionality.


Related commands
""""""""""""""""

:doc:`fix smd <fix_smd>`, :doc:`fix spring <fix_spring>`,
:doc:`fix plumed <fix_plumed>`

----------

.. _Fiorin:
**(Fiorin)** Fiorin, Klein, Henin, Mol. Phys. 111, 3345 (2013) https://doi.org/10.1080/00268976.2013.813594

.. _Colvars_LAMMPS_doc:
https://colvars.github.io/colvars-refman-lammps/colvars-refman-lammps.html
