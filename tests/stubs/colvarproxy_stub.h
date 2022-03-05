// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_STANDALONE_H
#define COLVARPROXY_STANDALONE_H

#define COLVARPROXY_VERSION COLVARS_VERSION


// Non-functional implementation that doesn't raise errors when critical
// functions aren't implemented
class colvarproxy_stub : public colvarproxy {

public:

  colvarproxy_stub();

  ~colvarproxy_stub() override;

  int setup() override;

  void log(std::string const &message) override;

  void error(std::string const &message) override;

  int set_unit_system(std::string const &units_in, bool check_only) override;

  int init_atom(int atom_number) override;

  int check_atom_id(int atom_number) override;

};


#endif
