cmake_minimum_required(VERSION 2.8.12)

project(colvars CXX)
get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

#set(CMAKE_VERBOSE_MAKEFILE on)

option(ENABLE_LEPTON "Build and link the Lepton library" OFF)

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/src/[^.]*.cpp)

if(ENABLE_LEPTON)
  set(LEPTON_DIR ${COLVARS_SOURCE_DIR}/lepton)
  set(CMAKE_CXX_STANDARD 11)
  file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
  add_library(lepton STATIC ${LEPTON_SOURCES})
  target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
endif()

add_library(colvars STATIC ${COLVARS_SOURCES})
target_compile_options(colvars PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-pedantic -Wall>)

if(ENABLE_LEPTON)
  target_compile_options(colvars PRIVATE -DLEPTON)
  target_include_directories(colvars PRIVATE ${LEPTON_DIR}/include)
endif()

message(STATUS "NOTE: this CMake recipe is useful for compilation tests only.")
message(STATUS "      To use the library, please follow instead the build recipes of the packages that include it.")
