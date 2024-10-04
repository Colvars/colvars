cmake_minimum_required(VERSION 3.13)

project(colvars CXX)
get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_VERBOSE_MAKEFILE ON)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

option(WARNINGS_ARE_ERRORS "Treats warnings as errors" OFF)

option(COLVARS_DEBUG "Build debugging binaries" OFF)

set(CMAKE_ADDRESS_SANITIZER OFF)
if(DEFINED ENV{CMAKE_ADDRESS_SANITIZER})
  set(CMAKE_ADDRESS_SANITIZER $ENV{CMAKE_ADDRESS_SANITIZER})
endif()
option(ADDRESS_SANITIZER "Build with Address Sanitizer for memory debugging" ${CMAKE_ADDRESS_SANITIZER})
if(ADDRESS_SANITIZER)
  add_compile_options(-fsanitize=address)
  add_link_options(-fsanitize=address)
endif()

if(NOT DEFINED LEPTON_DIR)
  set(LEPTON_DIR ${COLVARS_SOURCE_DIR}/lepton)
endif()

if(EXISTS ${LEPTON_DIR})
  set(COLVARS_LEPTON_DEFAULT ON)
else()
  set(COLVARS_LEPTON_DEFAULT OFF)
endif()

option(COLVARS_LEPTON "Build and link the Lepton library" ${COLVARS_LEPTON_DEFAULT})

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/src/[^.]*.cpp)

add_library(colvars ${COLVARS_SOURCES})
if(BUILD_SHARED_LIBS)
  set_property(TARGET colvars PROPERTY POSITION_INDEPENDENT_CODE 1)
endif()

target_compile_options(colvars PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-Wall -pedantic>)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION GREATER 6)
    target_compile_options(colvars PRIVATE -Wshadow-field)
  endif()
endif()

# # Note: disabled +w flags for Sun compiler due to their inconsistency
# target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:SunPro>:+w>)
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


if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set(COLVARS_OPENMP_DEFAULT OFF)
else()
  set(COLVARS_OPENMP_DEFAULT ON)
endif()
option(COLVARS_OPENMP "Build Colvars with (still limited) OpenMP support" ${COLVARS_OPENMP_DEFAULT})

if(COLVARS_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(colvars PRIVATE OpenMP::OpenMP_CXX)
  endif()
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
  target_link_libraries(colvars PUBLIC ${TCL_LIBRARY})
endif()

option(COLVARS_TORCH "Link against the Torch library" OFF)

if(COLVARS_TORCH)
  if (CMAKE_CXX_STANDARD LESS 17)
    message(STATUS "Increase the CMAKE_CXX_STANDARD to 17")
    set_property(TARGET colvars PROPERTY CXX_STANDARD 17)
  else()
    message(STATUS "Keep the CMAKE_CXX_STANDARD as ${CMAKE_CXX_STANDARD}")
  endif()
  if(DEFINED LIBTORCH_PREFIX)
    find_package(Torch REQUIRED PATHS file(JOIN ${LIBTORCH_PREFIX} "share/cmake/Torch"))
  else()
    find_package(Torch REQUIRED)
  endif()
  if(Torch_FOUND)
    target_compile_definitions(colvars PRIVATE -DCOLVARS_TORCH)
    target_compile_options(colvars PRIVATE ${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS})
    target_include_directories(colvars PRIVATE ${TORCH_INCLUDE_DIRS})
    target_link_libraries(colvars PUBLIC ${TORCH_LIBRARIES})
  endif()
endif()

option(BUILD_TOOLS "Build standalone tools" ON)

option(BUILD_TESTS "Build tests" ON)
 
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # TODO fix linkage errors with Visual Studio...
  set(BUILD_TOOLS OFF)
  set(BUILD_TESTS OFF)
endif()

option(BUILD_UNITTESTS "Build unit tests" ${BUILD_TESTS})

if(BUILD_TOOLS)
  add_subdirectory(${COLVARS_SOURCE_DIR}/colvartools colvartools)
endif()

if(COLVARS_DEBUG)
  target_compile_options(colvars PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-DCOLVARS_DEBUG -g -O0>)
  target_compile_options(colvars PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/DCOLVARS_DEBUG /DEBUG:FULL>)
endif()

set(CMAKE_COVERAGE OFF)
if(DEFINED ENV{CMAKE_COVERAGE})
  set(CMAKE_COVERAGE $ENV{CMAKE_COVERAGE})
endif()
option(ENABLE_COVERAGE "Enable code coverage reporting" ${CMAKE_COVERAGE})
if(ENABLE_COVERAGE AND (NOT CMAKE_CXX_COMPILER_ID STREQUAL "SunPro") AND (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC"))
  target_compile_options(colvars PUBLIC --coverage)
  target_link_options(colvars PUBLIC --coverage)
endif()


include(buildColvarsTests)
