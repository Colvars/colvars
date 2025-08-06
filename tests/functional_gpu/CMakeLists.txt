cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
include(CheckIncludeFile)
project(run_colvars_test_cuda LANGUAGES CXX)

if(NOT EXISTS ${COLVARS_SOURCE_DIR})
  set(COLVARS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../src" CACHE STRING "Colvars source code directory")
endif()

if(NOT DEFINED CMAKE_CUDA_STANDARD)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES native)
endif()

find_package(CUDAToolkit)
include(CheckLanguage)
check_language(CUDA)
enable_language(CUDA)

set(COLVARS_CUDA_DIR ${COLVARS_SOURCE_DIR}/cuda)

set(COLVARS_CUDA_HEADER ${COLVARS_CUDA_DIR}/colvaratoms_kernel.h
                        ${COLVARS_CUDA_DIR}/colvar_rotation_derivative_kernel.h
                        ${COLVARS_CUDA_DIR}/colvartypes_kernel.h)
set(COLVARS_CUDA_SOURCE ${COLVARS_CUDA_DIR}/colvaratoms_kernel.cu
                        ${COLVARS_CUDA_DIR}/colvar_rotation_derivative_kernel.cu
                        ${COLVARS_CUDA_DIR}/colvartypes_kernel.cu)

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/[^.]*.cpp)

add_definitions(-DCOLVARS_CUDA -DCOLVARS_DEBUG)

add_executable(
  run_colvars_test_cuda
  run_colvars_test_cuda.cpp
  ${COLVARS_SOURCES}
  ${COLVARS_CUDA_HEADER}
  ${COLVARS_CUDA_SOURCE})
target_include_directories(run_colvars_test_cuda
  PUBLIC "${COLVARS_SOURCE_DIR}"
  PUBLIC "${COLVARS_CUDA_DIR}"
  PUBLIC "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

if(CUDAToolkit_VERSION_MAJOR GREATER_EQUAL 10)
  target_link_libraries(run_colvars_test_cuda CUDA::cudart CUDA::nvtx3)
else()
  target_link_libraries(run_colvars_test_cuda CUDA::cudart CUDA::nvToolsExt)
endif()
set_property(TARGET run_colvars_test_cuda PROPERTY LANGUAGE CUDA)
# set_property(TARGET run_colvars_test_cuda PROPERTY CUDA_ARCHITECTURES OFF)
