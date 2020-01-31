// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


CVSCRIPT(bias_energy,
         "Get the current energy of this bias",
         0, 0,
         {},
         script->set_result_str(cvm::to_str(this_bias->get_energy()));
         return COLVARS_OK;
         )
