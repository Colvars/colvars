cmake_minimum_required(VERSION 3.1.3)

project(colvars CXX)
get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_VERBOSE_MAKEFILE on)

if(NOT DEFINED LEPTON_DIR)
  set(LEPTON_DIR ${COLVARS_SOURCE_DIR}/lepton)
else()
  if(NOT EXISTS ${LEPTON_DIR})
    message(FATAL_ERROR "Could not access ${LEPTON_DIR}")
  endif()
endif()

if(EXISTS ${LEPTON_DIR})
  set(COLVARS_LEPTON_DEFAULT ON)
else()
  set(COLVARS_LEPTON_DEFAULT OFF)
endif()

option(COLVARS_LEPTON "Build and link the Lepton library" ${COLVARS_LEPTON_DEFAULT})

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/src/[^.]*.cpp)

if(COLVARS_LEPTON)
  if(NOT EXISTS ${LEPTON_DIR})
    message(FATAL_ERROR "With -DCOLVARS_LEPTON=ON, the lepton folder must be copied into ${COLVARS_SOURCE_DIR} or provided by -DLEPTON_DIR=xxx")
  endif()
  set(CMAKE_CXX_STANDARD 11)
  file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
  add_library(lepton SHARED ${LEPTON_SOURCES})
  target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
endif()

add_library(colvars SHARED ${COLVARS_SOURCES})
add_library(colvars_static STATIC ${COLVARS_SOURCES})
set_property(TARGET colvars_static PROPERTY OUTPUT_NAME colvars)

target_compile_options(colvars PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-Wall -pedantic>)
target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-long-long>)
target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-c++11-long-long>)

if(COLVARS_LEPTON)
  target_compile_options(colvars PRIVATE -DLEPTON)
  target_include_directories(colvars PRIVATE ${LEPTON_DIR}/include)
  target_link_libraries(colvars lepton)
  # Silence warning
  target_compile_options(lepton PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-tautological-undefined-compare -Wno-unknown-warning-option>)
endif()

option(COLVARS_TCL "Link against the Tcl library" OFF)

if(COLVARS_TCL)
  find_library(TCL_LIBRARY NAMES tcl PATHS HINTS "${TCL_DIR}" PATH_SUFFIXES lib lib64)
  target_include_directories(colvars PRIVATE ${TCL_DIR}/include)
  target_compile_options(colvars PRIVATE -DCOLVARS_TCL)
  target_link_libraries(colvars tcl)
endif()

option(BUILD_TOOLS "Build standalone tools" ON)
if(BUILD_TOOLS)
  add_subdirectory(${COLVARS_SOURCE_DIR}/colvartools ./colvartools)
endif()

option(BUILD_UNITTESTS "Build unit tests" ON)
if(BUILD_UNITTESTS)
  add_subdirectory(${COLVARS_SOURCE_DIR}/tests/unittests ./unittests)
endif()

message(STATUS "NOTE: this CMake recipe is useful for compilation tests only.")
message(STATUS "      To use the library, please follow instead the build recipes of the packages that include it.")
