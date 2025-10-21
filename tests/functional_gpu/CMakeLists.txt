cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
include(CheckIncludeFile)
project(run_colvars_test_cuda LANGUAGES CXX)

option(BUILD_TESTS "Build tests" ON)
option(COLVARS_LEPTON "Build Colvars with Lepton" ON)

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
                        ${COLVARS_CUDA_DIR}/colvarcomp_distance_kernel.h
                        ${COLVARS_CUDA_DIR}/colvartypes_kernel.h)
set(COLVARS_CUDA_SOURCE ${COLVARS_CUDA_DIR}/colvaratoms_kernel.cu
                        ${COLVARS_CUDA_DIR}/colvar_rotation_derivative_kernel.cu
                        ${COLVARS_CUDA_DIR}/colvarcomp_distance_kernel.cu
                        ${COLVARS_CUDA_DIR}/colvartypes_kernel.cu)

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/src/[^.]*.cpp)

# Build Lepton
if(COLVARS_LEPTON)
  if(NOT DEFINED LEPTON_DIR)
    set(LEPTON_DIR "${COLVARS_SOURCE_DIR}/openmm-source/libraries/lepton")
    if(NOT EXISTS ${LEPTON_DIR})
      # Try the parent folder
      get_filename_component(LEPTON_DIR ${COLVARS_SOURCE_DIR} DIRECTORY)
      set(LEPTON_DIR "${LEPTON_DIR}/openmm-source/libraries/lepton")
    endif()
    if(NOT EXISTS ${LEPTON_DIR})
      # Giving up, cloning OpenMM into a sub-folder
      execute_process(COMMAND git clone --depth=1 https://github.com/openmm/openmm.git "${COLVARS_SOURCE_DIR}/openmm-source")
      set(LEPTON_DIR "${COLVARS_SOURCE_DIR}/openmm-source/libraries/lepton")
    endif()
    message(STATUS "Using Lepton library from: ${LEPTON_DIR}")
  endif()
  file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
  add_library(lepton OBJECT ${LEPTON_SOURCES})
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(lepton PRIVATE /DLEPTON_BUILDING_STATIC_LIBRARY /wd4267)
  endif()
  set_property(TARGET lepton PROPERTY POSITION_INDEPENDENT_CODE 1)
  target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
  # Silence warnings for Lepton
  target_compile_options(lepton PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-tautological-undefined-compare -Wno-unknown-warning-option>)
endif()

add_definitions(-DCOLVARS_CUDA)
if(COLVARS_LEPTON)
  add_definitions(-DLEPTON)
endif()

add_executable(
  run_colvars_test_cuda
  CLI11.hpp
  run_colvars_test_cuda.cpp
  ${COLVARS_SOURCES}
  ${COLVARS_CUDA_HEADER}
  ${COLVARS_CUDA_SOURCE})
target_include_directories(run_colvars_test_cuda
  PUBLIC "${COLVARS_SOURCE_DIR}/src"
  PUBLIC "${COLVARS_CUDA_DIR}"
  PUBLIC "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
  PUBLIC $<$<BOOL:${COLVARS_LEPTON}>:${LEPTON_DIR}/include>)

if(CUDAToolkit_VERSION_MAJOR GREATER_EQUAL 10)
  target_link_libraries(run_colvars_test_cuda PUBLIC CUDA::cudart PUBLIC CUDA::nvtx3)
else()
  target_link_libraries(run_colvars_test_cuda PUBLIC CUDA::cudart PUBLIC CUDA::nvToolsExt)
endif()

if(COLVARS_LEPTON)
  target_link_libraries(run_colvars_test_cuda PUBLIC lepton)
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
        COMMAND run_colvars_test_cuda --configuration_file ${TEST_NAME}/test.in --trajectory_file trajectory.xyz --output_prefix "${TEST_NAME}/test_out" --force
        WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/tests/functional_gpu)
      add_test(NAME "${TEST_NAME}_spiff"
        COMMAND bash compare_test.sh ${TEST_NAME}
        WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/tests/functional_gpu)
    endforeach()

    # Copy input files (coordinates, index files, etc)
    add_custom_command(TARGET copy_test_files POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${COLVARS_SOURCE_DIR}/tests/input_files
        ${CMAKE_BINARY_DIR}/tests/functional_gpu
      COMMENT "Copying test input files"
    )

    # Make the main target depend on the copy target
    add_dependencies(run_colvars_test_cuda copy_test_files)
  endif()
endif()

if(BUILD_TESTS)
  enable_testing()
endif()
