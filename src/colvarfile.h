// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef COLVARFILE_H
#define COLVARFILE_H

// OpenMM's Python binding supplies UTF-8 filenames. Use wide Windows I/O
// at the filesystem boundary. Other engines keep their native byte paths.
#ifdef COLVARS_UTF8_WINDOWS_PATHS
#include <filesystem>
inline std::filesystem::path colvars_file_path(char const *name)
{
  return std::filesystem::u8path(name);
}
#else
inline char const *colvars_file_path(char const *name)
{
  return name;
}
#endif

#endif
