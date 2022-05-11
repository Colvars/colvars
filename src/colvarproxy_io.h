// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_IO_H
#define COLVARPROXY_IO_H

#include <string>
#include <iostream>


/// Methods for data input/output
class colvarproxy_io {

public:

  /// Constructor
  colvarproxy_io();

  /// Destructor
  virtual ~colvarproxy_io();

  /// \brief Save the current frame number in the argument given
  // Returns error code
  virtual int get_frame(long int &);

  /// \brief Set the current frame number (as well as colvarmodule::it)
  // Returns error code
  virtual int set_frame(long int);

  /// \brief Rename the given file, before overwriting it
  virtual int backup_file(char const *filename);

  /// \brief Rename the given file, before overwriting it
  inline int backup_file(std::string const &filename)
  {
    return backup_file(filename.c_str());
  }

  /// Remove the given file (on Windows only, rename to filename.old)
  virtual int remove_file(char const *filename);

  /// Remove the given file (on Windows only, rename to filename.old)
  inline int remove_file(std::string const &filename)
  {
    return remove_file(filename.c_str());
  }

  /// Rename the given file
  virtual int rename_file(char const *filename, char const *newfilename);

  /// Rename the given file
  inline int rename_file(std::string const &filename,
                         std::string const &newfilename)
  {
    return rename_file(filename.c_str(), newfilename.c_str());
  }

  /// Prefix of the input state file to be read next
  inline std::string & input_prefix()
  {
    return input_prefix_str;
  }

  /// Default prefix to be used for all output files (final configuration)
  inline std::string & output_prefix()
  {
    return output_prefix_str;
  }

  /// Prefix of the restart (checkpoint) file to be written next
  inline std::string & restart_output_prefix()
  {
    return restart_output_prefix_str;
  }

  /// Default restart frequency (as set by the simulation engine)
  inline int default_restart_frequency() const
  {
    return restart_frequency_engine;
  }

  /// Buffer from which the input state information may be read
  inline char const * & input_buffer()
  {
    return input_buffer_;
  }

protected:

  /// Prefix of the input state file to be read next
  std::string input_prefix_str;

  /// Default prefix to be used for all output files (final configuration)
  std::string output_prefix_str;

  /// Prefix of the restart (checkpoint) file to be written next
  std::string restart_output_prefix_str;

  /// How often the simulation engine will write its own restart
  int restart_frequency_engine;

  /// \brief Currently opened output files: by default, these are ofstream objects.
  /// Allows redefinition to implement different output mechanisms
  std::list<std::ostream *> output_files;
  /// \brief Identifiers for output_stream objects: by default, these are the names of the files
  std::list<std::string>    output_stream_names;

  /// Buffer from which the input state information may be read
  char const *input_buffer_;
};


#endif
