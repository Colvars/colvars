cmake_minimum_required(VERSION 3.23)

project(OpenMMColvars
    VERSION 0.1.0
    LANGUAGES CXX)

option(OPENMM_COLVARS_BUILD_PYTHON "Build Python bindings for OpenMM-Colvars" ON)
option(OPENMM_COLVARS_ENABLE_CUDA_KERNELS
    "Build CUDA-platform host/device kernels for AlchemicalPMEForce when OpenMM CUDA internal headers are available"
    OFF)
option(OPENMM_COLVARS_FORCE_DISABLE_TABF_ALPHA
    "Build without TABF proxy hooks even if the Colvars source tree provides them"
    OFF)
set(OPENMM_COLVARS_PYTHON_INSTALL_DIR "" CACHE PATH
    "Python site-packages directory for installing openmm_colvars; empty means use the active Python environment")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(OPENMM_SOURCE_DIR "" CACHE PATH "Optional OpenMM source tree used for development builds and bundled Lepton")
if (NOT OPENMM_SOURCE_DIR AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests_dependence_src/openmm/CMakeLists.txt")
    set(OPENMM_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests_dependence_src/openmm" CACHE PATH "" FORCE)
endif()

find_package(OpenMM QUIET)
if (OpenMM_FOUND AND NOT TARGET OpenMM::OpenMM)
    if (TARGET OpenMM)
        add_library(OpenMM::OpenMM ALIAS OpenMM)
    else()
        message(FATAL_ERROR "find_package(OpenMM) succeeded but did not provide an OpenMM target")
    endif()
endif()

if (NOT OpenMM_FOUND AND NOT TARGET OpenMM::OpenMM)
    find_path(OPENMM_INCLUDE_DIR openmm/Context.h
        HINTS
            ${CMAKE_PREFIX_PATH}
            "$ENV{CONDA_PREFIX}"
        PATH_SUFFIXES
            include
            Library/include)
    find_library(OPENMM_LIBRARY NAMES OpenMM
        HINTS
            ${CMAKE_PREFIX_PATH}
            "$ENV{CONDA_PREFIX}"
        PATH_SUFFIXES
            lib
            Library/lib)
    if (OPENMM_INCLUDE_DIR AND OPENMM_LIBRARY)
        add_library(OpenMM::OpenMM UNKNOWN IMPORTED)
        set_target_properties(OpenMM::OpenMM PROPERTIES
            IMPORTED_LOCATION "${OPENMM_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENMM_INCLUDE_DIR}")
    endif()
endif()

set(OPENMM_COLVARS_USING_OPENMM_SOURCE_BUILD OFF)
if (NOT TARGET OpenMM::OpenMM AND OPENMM_SOURCE_DIR AND EXISTS "${OPENMM_SOURCE_DIR}/CMakeLists.txt")
    set(OPENMM_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_STATIC_LIB ON CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_CUDA_LIB OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_OPENCL_LIB OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_HIP_LIB OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_CPU_LIB OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_AMOEBA_PLUGIN OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_RPMD_PLUGIN OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_DRUDE_PLUGIN OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_PME_PLUGIN OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "" FORCE)
    set(OPENMM_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    list(APPEND CMAKE_MODULE_PATH "${OPENMM_SOURCE_DIR}/cmake_modules")
    add_subdirectory("${OPENMM_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/openmm-src" EXCLUDE_FROM_ALL)
    file(GLOB OPENMM_SERIALIZATION_SOURCES CONFIGURE_DEPENDS "${OPENMM_SOURCE_DIR}/serialization/src/*.cpp")
    set_source_files_properties(${OPENMM_SERIALIZATION_SOURCES}
        DIRECTORY "${OPENMM_SOURCE_DIR}"
        PROPERTIES COMPILE_DEFINITIONS "OPENMM_BUILDING_SHARED_LIBRARY;IEEE_8087")
    add_library(OpenMM::OpenMM ALIAS OpenMM_static)
    set(OPENMM_COLVARS_USING_OPENMM_SOURCE_BUILD ON)
endif()

if (NOT TARGET OpenMM::OpenMM)
    message(FATAL_ERROR
        "Could not find OpenMM. Set OpenMM_DIR to a directory containing "
        "OpenMMConfig.cmake, set OPENMM_INCLUDE_DIR and OPENMM_LIBRARY, "
        "or set OPENMM_SOURCE_DIR to an OpenMM source tree.")
endif()

set(OPENMM_COLVARS_OPENMM_INCLUDE_HINTS)
if (OPENMM_INCLUDE_DIR)
    list(APPEND OPENMM_COLVARS_OPENMM_INCLUDE_HINTS "${OPENMM_INCLUDE_DIR}")
endif()
get_target_property(OPENMM_COLVARS_OPENMM_INTERFACE_INCLUDES
    OpenMM::OpenMM INTERFACE_INCLUDE_DIRECTORIES)
if (OPENMM_COLVARS_OPENMM_INTERFACE_INCLUDES)
    list(APPEND OPENMM_COLVARS_OPENMM_INCLUDE_HINTS
        ${OPENMM_COLVARS_OPENMM_INTERFACE_INCLUDES})
endif()
list(REMOVE_DUPLICATES OPENMM_COLVARS_OPENMM_INCLUDE_HINTS)

set(COLVARS_SOURCE_DIR "" CACHE PATH "Path to the Colvars source directory containing src/colvarmodule.cpp")
if (NOT COLVARS_SOURCE_DIR AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests_dependence_src/colvars/src/colvarmodule.cpp")
    set(COLVARS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests_dependence_src/colvars" CACHE PATH "" FORCE)
endif()
if (NOT COLVARS_SOURCE_DIR)
    message(FATAL_ERROR "Set COLVARS_SOURCE_DIR to the Colvars repository root, for example -DCOLVARS_SOURCE_DIR=../colvars")
endif()

set(COLVARS_CORE_DIR "${COLVARS_SOURCE_DIR}/src")
if (NOT EXISTS "${COLVARS_CORE_DIR}/colvarmodule.cpp")
    message(FATAL_ERROR "COLVARS_SOURCE_DIR must contain src/colvarmodule.cpp")
endif()

set(OPENMM_COLVARS_HAS_TABF_ALPHA OFF)
set(OPENMM_COLVARS_HAS_TABF_ALPHA_BOUNDARIES OFF)
if (OPENMM_COLVARS_FORCE_DISABLE_TABF_ALPHA)
    message(STATUS "Colvars TABF alpha backend hooks disabled by OPENMM_COLVARS_FORCE_DISABLE_TABF_ALPHA")
elseif (EXISTS "${COLVARS_CORE_DIR}/colvarproxy_system.h")
    file(READ "${COLVARS_CORE_DIR}/colvarproxy_system.h" OPENMM_COLVARS_PROXY_SYSTEM_H)
    string(FIND "${OPENMM_COLVARS_PROXY_SYSTEM_H}" "get_tabf_alpha" OPENMM_COLVARS_TABF_ALPHA_POS)
    if (NOT OPENMM_COLVARS_TABF_ALPHA_POS EQUAL -1)
        set(OPENMM_COLVARS_HAS_TABF_ALPHA ON)
    endif()
    string(FIND "${OPENMM_COLVARS_PROXY_SYSTEM_H}" "get_tabf_alpha_boundaries"
        OPENMM_COLVARS_TABF_ALPHA_BOUNDARIES_POS)
    if (NOT OPENMM_COLVARS_TABF_ALPHA_BOUNDARIES_POS EQUAL -1)
        set(OPENMM_COLVARS_HAS_TABF_ALPHA_BOUNDARIES ON)
    endif()
endif()
if (OPENMM_COLVARS_HAS_TABF_ALPHA)
    message(STATUS "Colvars TABF alpha backend hooks detected")
else()
    message(STATUS "Colvars TABF alpha backend hooks not detected; building OpenMM-Colvars without TABF proxy overrides")
endif()

file(GLOB COLVARS_SOURCES CONFIGURE_DEPENDS "${COLVARS_CORE_DIR}/[^.]*.cpp")

set(LEPTON_DIR "" CACHE PATH "Optional Lepton source tree containing include/Lepton.h and src/*.cpp")
if (NOT LEPTON_DIR AND OPENMM_SOURCE_DIR AND EXISTS "${OPENMM_SOURCE_DIR}/libraries/lepton/include/Lepton.h")
    set(LEPTON_DIR "${OPENMM_SOURCE_DIR}/libraries/lepton" CACHE PATH "" FORCE)
endif()
if (NOT LEPTON_DIR)
    find_path(LEPTON_INCLUDE_DIR Lepton.h)
    find_library(LEPTON_LIBRARY NAMES Lepton lepton)
endif()
if (LEPTON_DIR AND EXISTS "${LEPTON_DIR}/include/Lepton.h")
    file(GLOB OPENMM_COLVARS_LEPTON_SOURCES CONFIGURE_DEPENDS "${LEPTON_DIR}/src/[^.]*.cpp")
    add_library(OpenMMColvarsLepton OBJECT ${OPENMM_COLVARS_LEPTON_SOURCES})
    set_target_properties(OpenMMColvarsLepton PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_include_directories(OpenMMColvarsLepton PUBLIC "${LEPTON_DIR}/include")
    target_compile_definitions(OpenMMColvarsLepton
        PUBLIC LEPTON_USE_STATIC_LIBRARIES
        PRIVATE LEPTON_BUILDING_STATIC_LIBRARY)
    set(OPENMM_COLVARS_LEPTON_OBJECTS $<TARGET_OBJECTS:OpenMMColvarsLepton>)
    set(OPENMM_COLVARS_LEPTON_INCLUDE_DIR "${LEPTON_DIR}/include")
elseif (LEPTON_INCLUDE_DIR)
    set(OPENMM_COLVARS_LEPTON_INCLUDE_DIR "${LEPTON_INCLUDE_DIR}")
    if (LEPTON_LIBRARY)
        set(OPENMM_COLVARS_LEPTON_LIBRARIES "${LEPTON_LIBRARY}")
    else()
        message(FATAL_ERROR
            "Found Lepton headers but not a Lepton library. Set LEPTON_DIR "
            "to build Lepton from source, or set LEPTON_LIBRARY.")
    endif()
else()
    message(FATAL_ERROR
        "Could not find Lepton for Colvars customFunction support. "
        "Set LEPTON_DIR, set OPENMM_SOURCE_DIR, or install Lepton headers and library.")
endif()

set(OPENMM_COLVARS_OPENMM_SOURCE_INCLUDES)
set(OPENMM_COLVARS_OPENMM_DEFINITIONS)
set(OPENMM_COLVARS_REFERENCE_PME_SOURCES)
set(OPENMM_COLVARS_REFERENCE_PME_INCLUDE_DIR)
if (OPENMM_SOURCE_DIR AND EXISTS "${OPENMM_SOURCE_DIR}/platforms/reference/include/ReferencePME.h")
    list(APPEND OPENMM_COLVARS_OPENMM_SOURCE_INCLUDES
        "${OPENMM_SOURCE_DIR}/platforms/reference/include"
        "${OPENMM_SOURCE_DIR}/libraries/pocketfft/include")
    list(APPEND OPENMM_COLVARS_REFERENCE_PME_SOURCES
        "${OPENMM_SOURCE_DIR}/platforms/reference/src/SimTKReference/ReferencePME.cpp")
else()
    find_path(OPENMM_COLVARS_REFERENCE_PME_INCLUDE_DIR ReferencePME.h
        HINTS ${OPENMM_COLVARS_OPENMM_INCLUDE_HINTS}
        PATH_SUFFIXES openmm/reference platforms/reference/include)
    if (OPENMM_COLVARS_REFERENCE_PME_INCLUDE_DIR)
        list(APPEND OPENMM_COLVARS_OPENMM_SOURCE_INCLUDES
            "${OPENMM_COLVARS_REFERENCE_PME_INCLUDE_DIR}")
    endif()
endif()
if (OPENMM_COLVARS_USING_OPENMM_SOURCE_BUILD AND
    OPENMM_SOURCE_DIR AND EXISTS "${OPENMM_SOURCE_DIR}/openmmapi/include/openmm/Force.h")
    list(APPEND OPENMM_COLVARS_OPENMM_SOURCE_INCLUDES
        "${OPENMM_SOURCE_DIR}/openmmapi/include"
        "${OPENMM_SOURCE_DIR}/openmmapi/src"
        "${OPENMM_SOURCE_DIR}/olla/include"
        "${OPENMM_SOURCE_DIR}/platforms/reference/include"
        "${OPENMM_SOURCE_DIR}/serialization/include"
        "${OPENMM_SOURCE_DIR}/libraries/lepton/include")
    set(OPENMM_COLVARS_OPENMM_DEFINITIONS OPENMM_USE_STATIC_LIBRARIES)
endif()
list(REMOVE_DUPLICATES OPENMM_COLVARS_OPENMM_SOURCE_INCLUDES)

set(OPENMM_COLVARS_CUDA_INCLUDES)
set(OPENMM_COLVARS_CUDA_DEFINITIONS)
if (OPENMM_COLVARS_ENABLE_CUDA_KERNELS)
    set(OPENMM_COLVARS_CUDA_INTERNAL_INCLUDES)
    if (OPENMM_SOURCE_DIR AND EXISTS "${OPENMM_SOURCE_DIR}/platforms/cuda/include/CudaContext.h")
        list(APPEND OPENMM_COLVARS_CUDA_INTERNAL_INCLUDES
            "${OPENMM_SOURCE_DIR}/platforms/cuda/include"
            "${OPENMM_SOURCE_DIR}/platforms/common/include")
    else()
        find_path(OPENMM_COLVARS_CUDA_CONTEXT_INCLUDE_DIR openmm/cuda/CudaContext.h
            HINTS ${OPENMM_COLVARS_OPENMM_INCLUDE_HINTS})
        find_path(OPENMM_COLVARS_COMMON_FFT_INCLUDE_DIR openmm/common/FFT3D.h
            HINTS ${OPENMM_COLVARS_OPENMM_INCLUDE_HINTS})
        if (OPENMM_COLVARS_CUDA_CONTEXT_INCLUDE_DIR AND OPENMM_COLVARS_COMMON_FFT_INCLUDE_DIR)
            list(APPEND OPENMM_COLVARS_CUDA_INTERNAL_INCLUDES
                "${OPENMM_COLVARS_CUDA_CONTEXT_INCLUDE_DIR}"
                "${OPENMM_COLVARS_COMMON_FFT_INCLUDE_DIR}")
        endif()
    endif()
    if (NOT OPENMM_COLVARS_CUDA_INTERNAL_INCLUDES)
        message(FATAL_ERROR
            "OPENMM_COLVARS_ENABLE_CUDA_KERNELS requires OpenMM CUDA internal "
            "headers. Provide OPENMM_SOURCE_DIR with platforms/cuda/include/"
            "CudaContext.h, or an OpenMM installation whose include directory "
            "contains openmm/cuda/CudaContext.h and openmm/common/FFT3D.h.")
    endif()
    set(OPENMM_CUDA_LIBRARY "" CACHE FILEPATH
        "OpenMMCUDA library required when OPENMM_COLVARS_ENABLE_CUDA_KERNELS is ON")
    if (NOT OPENMM_CUDA_LIBRARY)
        get_target_property(OPENMM_COLVARS_OPENMM_LOCATION OpenMM::OpenMM IMPORTED_LOCATION)
        if (OPENMM_COLVARS_OPENMM_LOCATION)
            get_filename_component(OPENMM_COLVARS_OPENMM_LIBRARY_DIR "${OPENMM_COLVARS_OPENMM_LOCATION}" DIRECTORY)
            find_library(OPENMM_CUDA_LIBRARY NAMES OpenMMCUDA
                HINTS
                    "${OPENMM_COLVARS_OPENMM_LIBRARY_DIR}"
                    "${OPENMM_COLVARS_OPENMM_LIBRARY_DIR}/plugins"
                NO_DEFAULT_PATH)
        endif()
        if (NOT OPENMM_CUDA_LIBRARY)
            find_library(OPENMM_CUDA_LIBRARY NAMES OpenMMCUDA
                HINTS
                    "$ENV{CONDA_PREFIX}/lib/plugins"
                    "$ENV{CONDA_PREFIX}/Library/lib/plugins")
        endif()
    endif()
    if (NOT OPENMM_CUDA_LIBRARY)
        message(FATAL_ERROR
            "OPENMM_COLVARS_ENABLE_CUDA_KERNELS requires the OpenMMCUDA library. "
            "Set OPENMM_CUDA_LIBRARY to libOpenMMCUDA.so/OpenMMCUDA.lib from the "
            "same OpenMM installation used by OpenMM::OpenMM.")
    endif()
    find_path(OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR cuda.h
        HINTS
            "$ENV{CUDA_HOME}/include"
            "$ENV{CUDA_PATH}/include"
            "$ENV{CONDA_PREFIX}/include")
    if (NOT OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR OR
        NOT EXISTS "${OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR}/builtin_types.h")
        message(FATAL_ERROR
            "OPENMM_COLVARS_ENABLE_CUDA_KERNELS requires CUDA driver headers "
            "(cuda.h and builtin_types.h). Set OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR "
            "to the CUDA include directory.")
    endif()
    list(APPEND OPENMM_COLVARS_CUDA_INCLUDES
        ${OPENMM_COLVARS_CUDA_INTERNAL_INCLUDES}
        "${OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR}")
    list(APPEND OPENMM_COLVARS_CUDA_DEFINITIONS OPENMM_COLVARS_ENABLE_CUDA_KERNELS=1)
endif()

if (NOT OPENMM_COLVARS_REFERENCE_PME_SOURCES AND
    NOT OPENMM_COLVARS_REFERENCE_PME_INCLUDE_DIR)
    message(FATAL_ERROR
        "AlchemicalPMEForce requires OpenMM ReferencePME. Provide an "
        "OPENMM_SOURCE_DIR with platforms/reference/include/ReferencePME.h, "
        "or an OpenMM installation whose include directory contains "
        "openmm/reference/ReferencePME.h.")
endif()

# Conditionally include CUDA kernel source
set(OPENMM_COLVARS_SOURCES
    src/alchemical_pme_force.cpp
    src/alchemical_pme_force_impl.cpp
    src/reference_alchemical_pme_kernel.cpp
    src/colvars_force.cpp
    src/colvars_force_impl.cpp
    src/colvars_tabf_integrator.cpp
    src/colvarproxy_tcl_stub.cpp
    src/colvarscript_stub.cpp
    src/openmm_colvar_proxy.cpp
    ${OPENMM_COLVARS_REFERENCE_PME_SOURCES}
    ${OPENMM_COLVARS_LEPTON_OBJECTS}
    ${COLVARS_SOURCES})

if(OPENMM_COLVARS_ENABLE_CUDA_KERNELS)
    list(APPEND OPENMM_COLVARS_SOURCES src/cuda_alchemical_pme_kernel.cpp)
endif()

add_library(OpenMMColvars STATIC ${OPENMM_COLVARS_SOURCES})

set_target_properties(OpenMMColvars PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(OpenMMColvars
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        ${OPENMM_COLVARS_OPENMM_SOURCE_INCLUDES}
        ${OPENMM_COLVARS_CUDA_INCLUDES}
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
        "${COLVARS_CORE_DIR}"
        "${OPENMM_COLVARS_LEPTON_INCLUDE_DIR}")

target_link_libraries(OpenMMColvars PUBLIC OpenMM::OpenMM)
target_link_libraries(OpenMMColvars PRIVATE ${OPENMM_COLVARS_LEPTON_LIBRARIES})
if (OPENMM_COLVARS_ENABLE_CUDA_KERNELS)
    target_link_libraries(OpenMMColvars PRIVATE "${OPENMM_CUDA_LIBRARY}")
endif()

target_compile_definitions(OpenMMColvars
    PUBLIC ${OPENMM_COLVARS_OPENMM_DEFINITIONS}
    PRIVATE COLVARS_DISABLE_TCL LEPTON ${OPENMM_COLVARS_CUDA_DEFINITIONS})
if (OPENMM_COLVARS_HAS_TABF_ALPHA)
    target_compile_definitions(OpenMMColvars PRIVATE OPENMM_COLVARS_HAS_TABF_ALPHA=1)
endif()
if (OPENMM_COLVARS_HAS_TABF_ALPHA_BOUNDARIES)
    target_compile_definitions(OpenMMColvars PRIVATE OPENMM_COLVARS_HAS_TABF_ALPHA_BOUNDARIES=1)
endif()
if (TARGET OpenMMColvarsLepton)
    target_link_libraries(OpenMMColvars PRIVATE OpenMMColvarsLepton)
endif()
if (MSVC)
    target_compile_options(OpenMMColvars PRIVATE /FS)
endif()

if (OPENMM_COLVARS_BUILD_PYTHON)
    find_package(Python COMPONENTS Interpreter Development.Module QUIET)
    find_program(SWIG_EXECUTABLE swig)

    if (Python_Interpreter_FOUND AND Python_Development.Module_FOUND AND SWIG_EXECUTABLE)
        set(OPENMM_COLVARS_PYTHON_DIR "${CMAKE_BINARY_DIR}/python")
        set(OPENMM_COLVARS_SWIG_WRAPPER "${OPENMM_COLVARS_PYTHON_DIR}/openmm_colvars_wrap.cpp")
        set(OPENMM_COLVARS_PYTHON_SHADOW "${OPENMM_COLVARS_PYTHON_DIR}/openmm_colvars.py")
        file(MAKE_DIRECTORY "${OPENMM_COLVARS_PYTHON_DIR}")

        if (NOT OPENMM_COLVARS_PYTHON_INSTALL_DIR)
            execute_process(
                COMMAND "${Python_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib'))"
                OUTPUT_VARIABLE OPENMM_COLVARS_DETECTED_PYTHON_INSTALL_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE
                COMMAND_ERROR_IS_FATAL ANY)
            set(OPENMM_COLVARS_PYTHON_INSTALL_DIR
                "${OPENMM_COLVARS_DETECTED_PYTHON_INSTALL_DIR}"
                CACHE PATH
                "Python site-packages directory for installing openmm_colvars"
                FORCE)
        endif()

        set(OPENMM_COLVARS_SWIG_INCLUDE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/python"
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
            "${COLVARS_CORE_DIR}"
            ${OPENMM_COLVARS_OPENMM_SOURCE_INCLUDES})
        if (OPENMM_INCLUDE_DIR)
            list(APPEND OPENMM_COLVARS_SWIG_INCLUDE_DIRS "${OPENMM_INCLUDE_DIR}")
        endif()
        get_target_property(OPENMM_COLVARS_OPENMM_TARGET_INCLUDES OpenMM::OpenMM INTERFACE_INCLUDE_DIRECTORIES)
        if (OPENMM_COLVARS_OPENMM_TARGET_INCLUDES)
            list(APPEND OPENMM_COLVARS_SWIG_INCLUDE_DIRS ${OPENMM_COLVARS_OPENMM_TARGET_INCLUDES})
        endif()
        list(REMOVE_DUPLICATES OPENMM_COLVARS_SWIG_INCLUDE_DIRS)
        set(OPENMM_COLVARS_SWIG_INCLUDE_ARGS)
        foreach(include_dir IN LISTS OPENMM_COLVARS_SWIG_INCLUDE_DIRS)
            list(APPEND OPENMM_COLVARS_SWIG_INCLUDE_ARGS "-I${include_dir}")
        endforeach()

        add_custom_command(
            OUTPUT
                "${OPENMM_COLVARS_SWIG_WRAPPER}"
                "${OPENMM_COLVARS_PYTHON_SHADOW}"
            COMMAND
                "${SWIG_EXECUTABLE}"
                -c++
                -python
                -py3
                -module openmm_colvars
                -outdir "${OPENMM_COLVARS_PYTHON_DIR}"
                -o "${OPENMM_COLVARS_SWIG_WRAPPER}"
                ${OPENMM_COLVARS_SWIG_INCLUDE_ARGS}
                "${CMAKE_CURRENT_SOURCE_DIR}/python/openmm_colvars.i"
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/python/openmm_colvars.i"
                "${CMAKE_CURRENT_SOURCE_DIR}/python/openmm_colvars_openmm_import.i"
                "${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/ColvarsForce.h"
                "${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/ColvarsTABFIntegrator.h"
                "${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/AlchemicalPMEForce.h"
            VERBATIM)

        Python_add_library(_openmm_colvars MODULE "${OPENMM_COLVARS_SWIG_WRAPPER}")
        set_target_properties(_openmm_colvars PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY "${OPENMM_COLVARS_PYTHON_DIR}"
            RUNTIME_OUTPUT_DIRECTORY "${OPENMM_COLVARS_PYTHON_DIR}"
            ARCHIVE_OUTPUT_DIRECTORY "${OPENMM_COLVARS_PYTHON_DIR}"
            OUTPUT_NAME "_openmm_colvars")
        target_include_directories(_openmm_colvars
            PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/include"
                ${OPENMM_COLVARS_OPENMM_SOURCE_INCLUDES}
                "${OPENMM_COLVARS_PYTHON_DIR}")
        target_link_libraries(_openmm_colvars PRIVATE OpenMMColvars Python::Module)
        target_compile_definitions(_openmm_colvars PRIVATE ${OPENMM_COLVARS_OPENMM_DEFINITIONS})

        add_custom_target(OpenMMColvarsPython ALL
            DEPENDS _openmm_colvars "${OPENMM_COLVARS_PYTHON_SHADOW}")

        install(TARGETS _openmm_colvars
            LIBRARY DESTINATION "${OPENMM_COLVARS_PYTHON_INSTALL_DIR}"
            RUNTIME DESTINATION "${OPENMM_COLVARS_PYTHON_INSTALL_DIR}")
        install(FILES "${OPENMM_COLVARS_PYTHON_SHADOW}"
            DESTINATION "${OPENMM_COLVARS_PYTHON_INSTALL_DIR}")
        install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/openmm_fep"
            DESTINATION "${OPENMM_COLVARS_PYTHON_INSTALL_DIR}"
            FILES_MATCHING PATTERN "*.py")
    else()
        message(WARNING
            "Python bindings were requested, but Python development headers or "
            "SWIG were not found. Configure with -DOPENMM_COLVARS_BUILD_PYTHON=OFF "
            "to silence this warning, or install Python development headers and SWIG.")
    endif()
endif()
