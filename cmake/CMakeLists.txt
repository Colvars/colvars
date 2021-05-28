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
  add_library(lepton STATIC ${LEPTON_SOURCES})
  set_property(TARGET lepton PROPERTY POSITION_INDEPENDENT_CODE 1)
  target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
endif()

# Object library to be used by both static and shared versions to avoid compiling every object twice
add_library(colvars_obj OBJECT ${COLVARS_SOURCES})
set_property(TARGET colvars_obj PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(colvars STATIC $<TARGET_OBJECTS:colvars_obj>)
add_library(colvars_shared SHARED $<TARGET_OBJECTS:colvars_obj>)
set_property(TARGET colvars_shared PROPERTY OUTPUT_NAME colvars)

target_compile_options(colvars_obj PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-Wall -pedantic>)
target_compile_options(colvars_obj PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-long-long>)
target_compile_options(colvars_obj PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-c++11-long-long>)

if(COLVARS_LEPTON)
  target_compile_options(colvars_obj PRIVATE -DLEPTON)
  target_include_directories(colvars_obj PRIVATE ${LEPTON_DIR}/include)
  target_link_libraries(colvars lepton)
  target_link_libraries(colvars_shared lepton)
  # Silence warning
  target_compile_options(lepton PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-tautological-undefined-compare -Wno-unknown-warning-option>)
endif()

option(COLVARS_TCL "Link against the Tcl library" OFF)

if(COLVARS_TCL)
  find_package(TCL)
  if(NOT ${TCL_FOUND})
    message(FATAL_ERROR "Tcl library not found.")
  endif()
  target_include_directories(colvars_obj PRIVATE ${TCL_INCLUDE_PATH})
  target_compile_options(colvars_obj PRIVATE -DCOLVARS_TCL)
  target_link_libraries(colvars tcl)
  target_link_libraries(colvars_shared tcl)
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
