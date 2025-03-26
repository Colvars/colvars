cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
include(CheckIncludeFile)
project(cudaglobalmastercolvars LANGUAGES CXX)

option(NAMD_TCL "If your NAMD is built with TCL, use ON" ON)
option(NAMD_NVTX_ENABLED "If your NAMD is built with -DNAMD_NVTX_ENABLED, use ON" OFF)
set(NAMD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE STRING "NAMD source code directory")
set(NAMD_HEADER_DIR "${NAMD_DIR}/src" CACHE STRING "NAMD header directory")
set(MOLFILE_PLUGIN_HEADER_DIR "${NAMD_DIR}/plugins/include" CACHE STRING "Molfile plugin header directory")
set(COLVARS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src" CACHE STRING "Colvars source code directory")
message(STATUS "NAMD_HEADER_DIR is ${NAMD_HEADER_DIR}")
message(STATUS "MOLFILE_PLUGIN_HEADER_DIR is ${MOLFILE_PLUGIN_HEADER_DIR}")
message(STATUS "COLVARS_SOURCE_DIR is ${COLVARS_SOURCE_DIR}")

if(NOT EXISTS ${NAMD_HEADER_DIR})
  message(FATAL_ERROR "The NAMD header directory ${NAMD_HEADER_DIR} does not exist!")
endif()

if(NOT EXISTS ${MOLFILE_PLUGIN_HEADER_DIR})
  message(FATAL_ERROR "The molfile plugin header directory ${MOLFILE_PLUGIN_HEADER_DIR} does not exist!")
endif()

if(NOT EXISTS ${COLVARS_SOURCE_DIR})
  message(FATAL_ERROR "The Colvars source directory ${COLVARS_SOURCE_DIR} does not exist!")
endif()

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/[^.]*.cpp)

find_package(CUDAToolkit)
include(CheckLanguage)
check_language(CUDA)
enable_language(CUDA)

if(NOT DEFINED CMAKE_CUDA_STANDARD)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

if(NAMD_TCL)
  message(STATUS "add -DNAMD_TCL to compiler flags")
  add_definitions(-DNAMD_TCL)
endif()

if(NAMD_NVTX_ENABLED)
  message(STATUS "add -DNAMD_NVTX_ENABLED and -DCUDAGLOBALMASTERCOLVARS_CUDA_PROFILING to compiler flags")
  add_definitions(-DNAMD_NVTX_ENABLED -DCUDAGLOBALMASTERCOLVARS_CUDA_PROFILING)
endif()

add_definitions(-DNAMD_CUDA -DNODEGROUP_FORCE_REGISTER -DCOLVARS_CUDA)
add_library(cudaglobalmastercolvars SHARED
            colvarproxy_cudaglobalmaster.h
            colvarproxy_cudaglobalmaster.C
            colvarproxy_cudaglobalmaster_kernel.h
            colvarproxy_cudaglobalmaster_kernel.cu
            ${COLVARS_SOURCES})
target_include_directories(cudaglobalmastercolvars PUBLIC "${COLVARS_SOURCE_DIR}" PUBLIC "${NAMD_HEADER_DIR}" PUBLIC "${MOLFILE_PLUGIN_HEADER_DIR}")
target_link_libraries(cudaglobalmastercolvars CUDA::cudart CUDA::nvToolsExt)
target_compile_options(cudaglobalmastercolvars PRIVATE -Wno-register)
set_property(TARGET cudaglobalmastercolvars PROPERTY LANGUAGE CUDA)
set_property(TARGET cudaglobalmastercolvars PROPERTY CUDA_ARCHITECTURES OFF)
