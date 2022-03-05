cmake_minimum_required(VERSION 3.1.3)

project(colvars CXX)
get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_VERBOSE_MAKEFILE ON)

option(WARNINGS_ARE_ERRORS "Treats warnings as errors" OFF)

if(NOT DEFINED LEPTON_DIR)
  set(LEPTON_DIR ${COLVARS_SOURCE_DIR}/lepton)
endif()

if(EXISTS ${LEPTON_DIR})
  set(COLVARS_LEPTON_DEFAULT ON)
else()
  set(COLVARS_LEPTON_DEFAULT OFF)
endif()

option(COLVARS_LEPTON "Build and link the Lepton library" ${COLVARS_LEPTON_DEFAULT})

if(COLVARS_LEPTON)
  if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  else()
    if(${CMAKE_CXX_STANDARD} STREQUAL "98")
      message(FATAL_ERROR "With -DCOLVARS_LEPTON=ON, CMAKE_CXX_STANDARD must be 11 or later")
    endif()
  endif()
endif()

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/src/[^.]*.cpp)

add_library(colvars ${COLVARS_SOURCES})
if(BUILD_SHARED_LIBS)
  set_property(TARGET colvars PROPERTY POSITION_INDEPENDENT_CODE 1)
endif()

target_compile_options(colvars PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-Wall -pedantic>)

target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-long-long>)
target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-c++11-long-long>)

target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:SunPro>:+w>)
# # Uncomment the line below to enable more thorough checks with Sun Compiler
# target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:SunPro>:+w2>)

# Tell MSVC to actually report the C++ version
target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>)
# Use "production level" warnings, but silence those known to be harmless
target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/W3 /wd4244 /wd4267>)

if(WARNINGS_ARE_ERRORS)
  target_compile_options(colvars PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-Werror>)
  target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:SunPro>:-xwe>)
  target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/WX>)
endif()

include(buildColvarsLepton)

option(COLVARS_TCL "Link against the Tcl library" OFF)

if(COLVARS_TCL)
  if(DEFINED TCL_DIR)
    message(STATUS "Using Tcl headers and libraries from: ${TCL_DIR}")
    if(NOT DEFINED TCL_LIBRARY)
      find_library(TCL_LIBRARY NAMES tcl PATHS HINTS "${TCL_DIR}" PATH_SUFFIXES lib lib64)
    else()
      message(STATUS "Using Tcl library from explicit path: ${TCL_LIBRARY}")
    endif()
    set(TCL_INCLUDE_PATH ${TCL_DIR}/include)
  else()
    find_package(TCL)
    if(NOT ${TCL_FOUND})
      message(FATAL_ERROR "Tcl package not found.")
    endif()
  endif()
  target_include_directories(colvars PRIVATE ${TCL_INCLUDE_PATH})
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(colvars PRIVATE /DCOLVARS_TCL)
  else()
    target_compile_options(colvars PRIVATE -DCOLVARS_TCL)
  endif()
  target_link_libraries(colvars ${TCL_LIBRARY})
endif()

option(BUILD_TOOLS "Build standalone tools" OFF)
option(BUILD_UNITTESTS "Build unit tests" OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # TODO fix linkage errors with Visual Studio...
  set(BUILD_TOOLS OFF)
  set(BUILD_UNITTESTS OFF)
endif()

if(BUILD_TOOLS)
  add_subdirectory(${COLVARS_SOURCE_DIR}/colvartools colvartools)
endif()

option(BUILD_TESTS "Build tests" ON)

include(buildColvarsTests)
