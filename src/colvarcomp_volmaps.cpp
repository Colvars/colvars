// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"



colvar::map_total::map_total()
  : cvc(), volmap_index(-1)
{
  function_type = "map_total";
  volmap_id = -1;
  volmap_index = -1;
  atoms = NULL;
  x.type(colvarvalue::type_scalar);
}


colvar::map_total::map_total(std::string const &conf)
  : cvc(), volmap_index(-1)
{
  function_type = "map_total";
  volmap_id = -1;
  volmap_index = -1;
  atoms = NULL;
  x.type(colvarvalue::type_scalar);
  map_total::init(conf);
}


int colvar::map_total::init(std::string const &conf)
{
  int error_code = cvc::init(conf);
  colvarproxy *proxy = cvm::main()->proxy;
  get_keyval(conf, "mapName", volmap_name, volmap_name);
  get_keyval(conf, "mapID", volmap_id, volmap_id);

  if ((volmap_name.size() > 0) && (volmap_id >= 0)) {
    error_code |=
      cvm::error("Error: mapName and mapID are mutually exclusive.\n");
  }

  // Parse optional group
  atoms = parse_group(conf, "atoms", true);
  if (atoms != NULL) {

    // Using internal selection
    if (volmap_name.size()) {
      error_code |= proxy->check_volmap_by_name(volmap_name);
    }
    if (volmap_id >= 0) {
      error_code |= proxy->check_volmap_by_id(volmap_id);
    }

  } else {

    // Using selection from the MD engine
    if (volmap_name.size()) {
      volmap_index = proxy->init_volmap_by_name(volmap_name);
    }
    if (volmap_id >= 0) {
      volmap_index = proxy->init_volmap_by_id(volmap_id);
    }
    error_code |= volmap_index > 0 ? COLVARS_OK : INPUT_ERROR;
  }

  if (volmap_name.size() > 0) {
    volmap_id = proxy->get_volmap_id_from_name(volmap_name.c_str());
  }

  return error_code;
}


void colvar::map_total::calc_value()
{
  colvarproxy const *proxy = cvm::main()->proxy;
  if (atoms != NULL) {
    x.real_value = 0.0;
    proxy->compute_volmap(volmap_id, atoms->begin(), atoms->end(),
                          &(x.real_value), NULL);
  } else {
    x.real_value = proxy->get_volmap_value(volmap_index);
  }
}


void colvar::map_total::calc_gradients()
{
  // Computed in calc_value() or by the MD engine
}


void colvar::map_total::apply_force(colvarvalue const &force)
{
  colvarproxy *proxy = cvm::main()->proxy;
  if (atoms) {
    if (!atoms->noforce)
      atoms->apply_colvar_force(force.real_value);
  } else {
    proxy->apply_volmap_force(volmap_index, force.real_value);
  }
}
