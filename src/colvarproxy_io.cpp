// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

// Using access() to check if a file exists (until we can assume C++14/17)
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#if defined(WIN32)
#include <io.h>
#endif

#include <cerrno>

#include <sstream>
#include <cstring>
#include <cstdio>

#include "colvarmodule.h"
#include "colvarproxy_io.h"


colvarproxy_io::colvarproxy_io()
{
  input_buffer_ = NULL;
  restart_frequency_engine = 0;
}


colvarproxy_io::~colvarproxy_io() {}


int colvarproxy_io::get_frame(long int&)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::set_frame(long int)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::backup_file(char const *filename)
{
  // Simplified version of NAMD_file_exists()
  int exit_code;
  do {
#if defined(WIN32) && !defined(__CYGWIN__)
    // We could use _access_s here, but it is probably too new
    exit_code = _access(filename, 00);
#else
    exit_code = access(filename, F_OK);
#endif
  } while ((exit_code != 0) && (errno == EINTR));
  if (exit_code != 0) {
    if (errno == ENOENT) {
      // File does not exist
      return COLVARS_OK;
    } else {
      return cvm::error("Unknown error while checking if file \""+
                        std::string(filename)+"\" exists.\n", COLVARS_ERROR);
    }
  }

  // The file exists, then rename it
  if (std::string(filename).rfind(std::string(".colvars.state")) !=
      std::string::npos) {
    return rename_file(filename, (std::string(filename)+".old").c_str());
  } else {
    return rename_file(filename, (std::string(filename)+".BAK").c_str());
  }
}


int colvarproxy_io::remove_file(char const *filename)
{
  int error_code = COLVARS_OK;
#if defined(WIN32) && !defined(__CYGWIN__)
  // Because the file may be open by other processes, rename it to filename.old
  std::string const renamed_file(std::string(filename)+".old");
  // It may still be there from an interrupted run, so remove it to be safe
  std::remove(renamed_file.c_str());
  int rename_exit_code = 0;
  while ((rename_exit_code = std::rename(filename,
                                         renamed_file.c_str())) != 0) {
    if (errno == EINTR) continue;
    error_code |= COLVARS_FILE_ERROR;
    break;
  }
  // Ask to remove filename.old, but ignore any errors raised
  std::remove(renamed_file.c_str());
#else
  if (std::remove(filename)) {
    if (errno != ENOENT) {
      error_code |= COLVARS_FILE_ERROR;
    }
  }
#endif
  if (error_code != COLVARS_OK) {
    return cvm::error("Error: in removing file \""+std::string(filename)+
                      "\".\n.",
                      error_code);
  }
  return COLVARS_OK;
}


int colvarproxy_io::rename_file(char const *filename, char const *newfilename)
{
  int error_code = COLVARS_OK;
#if defined(WIN32) && !defined(__CYGWIN__)
  // On straight Windows, must remove the destination before renaming it
  error_code |= remove_file(newfilename);
#endif
  int rename_exit_code = 0;
  while ((rename_exit_code = std::rename(filename, newfilename)) != 0) {
    if (errno == EINTR) continue;
    // Call log() instead of error to allow the next try
    cvm::log("Error: in renaming file \""+std::string(filename)+"\" to \""+
             std::string(newfilename)+"\".\n.");
    error_code |= COLVARS_FILE_ERROR;
    if (errno == EXDEV) continue;
    break;
  }
  return rename_exit_code ? error_code : COLVARS_OK;
}
