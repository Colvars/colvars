# General build script suitable for cross-platform CI tests

get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()

if(NOT DEFINED LEPTON_DIR)
  set(LEPTON_DIR "${COLVARS_SOURCE_DIR}/openmm-source/libraries/lepton")
  if(NOT EXISTS ${LEPTON_DIR})
    execute_process(COMMAND git clone --depth=1 https://github.com/openmm/openmm.git "${COLVARS_SOURCE_DIR}/openmm-source")
  endif()
  message(STATUS "Using Lepton library from: ${LEPTON_DIR}")
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND}
  -S cmake
  -B build
  -D CMAKE_BUILD_TYPE=RelWithDebinfo
  -D CMAKE_C_COMPILER_LAUNCHER=ccache
  -D CMAKE_CXX_COMPILER_LAUNCHER=ccache
  -D WARNINGS_ARE_ERRORS=ON
  -D CMAKE_VERBOSE_MAKEFILE=ON
  -D CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
  -D COLVARS_TCL=ON
  -D COLVARS_LEPTON=ON
  -D LEPTON_DIR=${LEPTON_DIR}
  RESULT_VARIABLE result
  )
if (NOT result EQUAL 0)
  message(FATAL_ERROR "Error generating CMake buildsystem.")
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND}
  --build build
  --parallel
  )
if (NOT result EQUAL 0)
  message(FATAL_ERROR "Error building library.")
endif()
