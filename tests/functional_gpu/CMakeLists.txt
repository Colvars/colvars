cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
include(CheckIncludeFile)
project(run_colvars_test_cuda LANGUAGES CXX)

option(BUILD_TESTS "Build tests" ON)
option(OUTPUT_FORCE "Output the applied forces for debugging" ON)

if(NOT EXISTS ${COLVARS_SOURCE_DIR})
  set(COLVARS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE STRING "Colvars source code directory")
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

set(COLVARS_CUDA_DIR ${COLVARS_SOURCE_DIR}/src/cuda)

set(COLVARS_CUDA_HEADER ${COLVARS_CUDA_DIR}/colvaratoms_kernel.h
                        ${COLVARS_CUDA_DIR}/colvar_rotation_derivative_kernel.h
                        ${COLVARS_CUDA_DIR}/colvartypes_kernel.h)
set(COLVARS_CUDA_SOURCE ${COLVARS_CUDA_DIR}/colvaratoms_kernel.cu
                        ${COLVARS_CUDA_DIR}/colvar_rotation_derivative_kernel.cu
                        ${COLVARS_CUDA_DIR}/colvartypes_kernel.cu)

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/src/[^.]*.cpp)

add_definitions(-DCOLVARS_CUDA)

if(OUTPUT_FORCE)
  add_definitions(-DSTUB_OUTPUT_FORCE)
endif()

add_executable(
  run_colvars_test_cuda
  run_colvars_test_cuda.cpp
  ${COLVARS_SOURCES}
  ${COLVARS_CUDA_HEADER}
  ${COLVARS_CUDA_SOURCE})
target_include_directories(run_colvars_test_cuda
  PUBLIC "${COLVARS_SOURCE_DIR}/src"
  PUBLIC "${COLVARS_CUDA_DIR}"
  PUBLIC "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

if(CUDAToolkit_VERSION_MAJOR GREATER_EQUAL 10)
  target_link_libraries(run_colvars_test_cuda CUDA::cudart CUDA::nvtx3)
else()
  target_link_libraries(run_colvars_test_cuda CUDA::cudart CUDA::nvToolsExt)
endif()
set_property(TARGET run_colvars_test_cuda PROPERTY LANGUAGE CUDA)
# set_property(TARGET run_colvars_test_cuda PROPERTY CUDA_ARCHITECTURES OFF)

if(BUILD_TESTS)
  if(UNIX)
    add_custom_target(copy_test_files ALL
      COMMENT "Copying test input files and configurations"
      VERBATIM
    )

    # Copy the Colvars configuration files
    file(GLOB TEST_CONFIG_FILES ${COLVARS_SOURCE_DIR}/tests/input_files/*/test.in)
    if(NOT COLVARS_TORCH)
      # TODO create a way to detect test dependencies at some point
      list(REMOVE_ITEM TEST_CONFIG_FILES ${COLVARS_SOURCE_DIR}/tests/input_files/torchann-dihedral_harmonic-fixed/test.in)
    endif()

    foreach(TEST_CONFIG_FILE ${TEST_CONFIG_FILES})
      get_filename_component(TEST_NAME ${TEST_CONFIG_FILE} DIRECTORY)
      get_filename_component(TEST_NAME ${TEST_NAME} NAME)
      add_test(NAME ${TEST_NAME}
        COMMAND run_colvars_test_cuda ${TEST_NAME}/test.in trajectory.xyz
        WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/tests/functional)
    endforeach()

    # Copy input files (coordinates, index files, etc)
    add_custom_command(TARGET copy_test_files POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${COLVARS_SOURCE_DIR}/tests/input_files
        ${CMAKE_BINARY_DIR}/tests/functional
      COMMENT "Copying test input files"
    )

    # Make the main target depend on the copy target
    add_dependencies(run_colvars_test_cuda copy_test_files)
  endif()
endif()

if(BUILD_TESTS)
  enable_testing()
endif()
