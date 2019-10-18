cmake_minimum_required(VERSION 3.1.3)

project(colvars CXX)
get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_VERBOSE_MAKEFILE on)

option(COLVARS_LEPTON "Build and link the Lepton library" OFF)

if(NOT DEFINED LEPTON_DIR)
  set(LEPTON_DIR ${COLVARS_SOURCE_DIR}/lepton)
else()
  if(NOT EXISTS ${LEPTON_DIR})
    message(FATAL_ERROR "Could not access ${LEPTON_DIR}")
  endif()
endif()
if(EXISTS ${LEPTON_DIR})
  set(COLVARS_LEPTON ON)
endif()

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/src/[^.]*.cpp)

if(COLVARS_LEPTON)
  if(NOT EXISTS ${LEPTON_DIR})
    message(FATAL_ERROR "With -DCOLVARS_LEPTON=ON, the lepton folder must be copied into ${COLVARS_SOURCE_DIR} or provided by -DLEPTON_DIR=xxx")
  endif()
  set(CMAKE_CXX_STANDARD 11)
  file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
  add_library(lepton STATIC ${LEPTON_SOURCES})
  target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
endif()

add_library(colvars STATIC ${COLVARS_SOURCES})

target_compile_options(colvars PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-Wall -pedantic>)
target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-long-long>)
target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-c++11-long-long>)

if(COLVARS_LEPTON)
  target_compile_options(colvars PRIVATE -DLEPTON)
  target_include_directories(colvars PRIVATE ${LEPTON_DIR}/include)
  # Silence warning
  target_compile_options(lepton PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-tautological-undefined-compare -Wno-unknown-warning-option>)
endif()

message(STATUS "NOTE: this CMake recipe is useful for compilation tests only.")
message(STATUS "      To use the library, please follow instead the build recipes of the packages that include it.")
