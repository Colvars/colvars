// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


CVSCRIPT(colvar_value,
         "Get the current value of this colvar",
         0, 0,
         "",
         script->set_result_str(this_colvar->value().to_simple_string());
         return COLVARS_OK;
         )
