# This file is part of the Collective Variables module (Colvars).
# The original version of Colvars and its updates are located at:
# https://github.com/Colvars/colvars
# Please update all Colvars source files before making any changes.
# If you wish to distribute your changes, please submit them to the
# Colvars repository at GitHub.

function(gmx_manage_colvars)
  file(GLOB COLVARS_SOURCES ${PROJECT_SOURCE_DIR}/src/external/colvars/*.cpp)
  add_library(colvars OBJECT ${COLVARS_SOURCES})
  # Colvars requires a correct definition of __cplusplus, which MSVC doesn't give by default
  target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>)
endfunction()

function(gmx_include_colvars_headers)
  target_include_directories(libgromacs PRIVATE ${PROJECT_SOURCE_DIR}/src/external/colvars)
endfunction()
