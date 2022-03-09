cmake_minimum_required(VERSION 3.1.3)

project(colvars CXX)
get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

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

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/src/[^.]*.cpp)


# Object library to be used by both static and shared versions to avoid compiling every object twice
add_library(colvars_obj OBJECT ${COLVARS_SOURCES})
set_property(TARGET colvars_obj PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(colvars STATIC $<TARGET_OBJECTS:colvars_obj>)
add_library(colvars_shared SHARED $<TARGET_OBJECTS:colvars_obj>)
set_property(TARGET colvars_shared PROPERTY OUTPUT_NAME colvars)

target_compile_options(colvars_obj PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-Wall -pedantic>)

target_compile_options(colvars_obj PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-long-long>)
target_compile_options(colvars_obj PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-c++11-long-long>)

target_compile_options(colvars_obj PRIVATE $<$<CXX_COMPILER_ID:SunPro>:+w>)
# # Uncomment the line below to enable more thorough checks with Sun Compiler
# target_compile_options(colvars_obj PRIVATE $<$<CXX_COMPILER_ID:SunPro>:+w2>)

# Tell MSVC to actually report the C++ version
target_compile_options(colvars_obj PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/W2 /Zc:__cplusplus>)

if(WARNINGS_ARE_ERRORS)
  target_compile_options(colvars_obj PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Intel>>:-Werror>)
  target_compile_options(colvars_obj PRIVATE $<$<CXX_COMPILER_ID:SunPro>:-xwe>)
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
  target_include_directories(colvars_obj PRIVATE ${TCL_INCLUDE_PATH})
  target_compile_options(colvars_obj PRIVATE -DCOLVARS_TCL)
  target_link_libraries(colvars ${TCL_LIBRARY})
  target_link_libraries(colvars_shared ${TCL_LIBRARY})
endif()

option(BUILD_TOOLS "Build standalone tools" ON)
if(BUILD_TOOLS)
  add_subdirectory(${COLVARS_SOURCE_DIR}/colvartools ./colvartools)
endif()

option(BUILD_UNITTESTS "Build unit tests" ON)
if(BUILD_UNITTESTS)
  add_subdirectory(${COLVARS_SOURCE_DIR}/tests/unittests ./unittests)
endif()

