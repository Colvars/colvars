cmake_minimum_required(VERSION 3.24)
project(OpenMMColvars LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(OPENMM_COLVARS_ENABLE_CUDA "Build the native CUDA/Colvars device bridge" OFF)
option(OPENMM_COLVARS_BUILD_PYTHON "Build Python bindings" ON)
set(OPENMM_ROOT "" CACHE PATH "OpenMM prefix containing include and lib")
set(OPENMM_SOURCE_DIR "" CACHE PATH "Matching OpenMM sources (for optional Lepton)")
set(OPENMM_COLVARS_PYTHON_INSTALL_DIR "" CACHE PATH "Python installation directory")
get_filename_component(COLVARS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(COLVARS_CORE_DIR "${COLVARS_SOURCE_DIR}/src")

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(COMMAND "${Python_EXECUTABLE}" -c
    "import pathlib,openmm; print(pathlib.Path(openmm.__file__).resolve().parent.parent/'OpenMM.libs')"
    RESULT_VARIABLE OPENMM_PYTHON_RESULT OUTPUT_VARIABLE OPENMM_PYTHON_ROOT
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
find_path(OPENMM_INCLUDE_DIR openmm/Context.h
    HINTS "${OPENMM_ROOT}" "${OPENMM_PYTHON_ROOT}" "$ENV{CONDA_PREFIX}"
    PATH_SUFFIXES include Library/include REQUIRED)
find_library(OPENMM_LIBRARY NAMES OpenMM
    HINTS "${OPENMM_ROOT}" "${OPENMM_PYTHON_ROOT}" "$ENV{CONDA_PREFIX}"
    PATH_SUFFIXES lib Library/lib REQUIRED)
get_filename_component(OPENMM_LIBRARY_DIR "${OPENMM_LIBRARY}" DIRECTORY)
message(STATUS "OpenMM headers: ${OPENMM_INCLUDE_DIR}")
message(STATUS "OpenMM library: ${OPENMM_LIBRARY}")
if(NOT EXISTS "${OPENMM_INCLUDE_DIR}/openmm/internal/CustomCPPForceImpl.h")
    message(FATAL_ERROR "OpenMM with CustomCPPForceImpl is required (OpenMM 8.1 or newer)")
endif()

file(GLOB COLVARS_CORE_SOURCES CONFIGURE_DEPENDS "${COLVARS_CORE_DIR}/*.cpp")
add_library(OpenMMColvars STATIC
    src/ColvarsForce.cpp src/ColvarsForceImpl.cpp src/colvarproxy_openmm.cpp
    ${COLVARS_CORE_SOURCES})
set_target_properties(OpenMMColvars PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(OpenMMColvars PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    PRIVATE "${COLVARS_CORE_DIR}")
target_include_directories(OpenMMColvars SYSTEM PUBLIC "${OPENMM_INCLUDE_DIR}")
target_link_libraries(OpenMMColvars PUBLIC "${OPENMM_LIBRARY}")
if(WIN32 AND NOT CYGWIN)
    target_compile_definitions(OpenMMColvars PRIVATE COLVARS_UTF8_WINDOWS_PATHS)
endif()
if(MSVC)
    target_compile_options(OpenMMColvars PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:/Zc:__cplusplus;/utf-8>")
else()
    target_compile_options(OpenMMColvars PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall>)
endif()

# Reuse the matching OpenMM source's Lepton implementation when supplied.
set(LEPTON_DIR "" CACHE PATH "Optional explicit Lepton source directory")
if(LEPTON_DIR)
    set(COLVARS_LEPTON_DIR "${LEPTON_DIR}")
else()
    # Keep the inferred path out of the cache: users may supply or change
    # OPENMM_SOURCE_DIR on a later configure of the same build directory.
    set(COLVARS_LEPTON_DIR "${OPENMM_SOURCE_DIR}/libraries/lepton")
endif()
if(EXISTS "${COLVARS_LEPTON_DIR}/include/Lepton.h")
    file(GLOB LEPTON_SOURCES CONFIGURE_DEPENDS "${COLVARS_LEPTON_DIR}/src/*.cpp")
    target_sources(OpenMMColvars PRIVATE ${LEPTON_SOURCES})
    target_include_directories(OpenMMColvars SYSTEM PRIVATE "${COLVARS_LEPTON_DIR}/include")
    target_compile_definitions(OpenMMColvars PRIVATE LEPTON LEPTON_BUILDING_STATIC_LIBRARY)
    message(STATUS "Colvars customFunction support enabled with ${COLVARS_LEPTON_DIR}")
else()
    message(STATUS "Colvars customFunction disabled: supply OPENMM_SOURCE_DIR or LEPTON_DIR to enable")
endif()

if(OPENMM_COLVARS_ENABLE_CUDA)
    find_package(CUDAToolkit 12.0 REQUIRED)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_library(OPENMM_CUDA_LIBRARY NAMES OpenMMCUDA
        HINTS "${OPENMM_LIBRARY_DIR}/plugins" "${OPENMM_LIBRARY_DIR}" REQUIRED)
    if(NOT EXISTS "${OPENMM_INCLUDE_DIR}/openmm/cuda/CudaContext.h")
        message(FATAL_ERROR "CUDA build requires installed OpenMM CUDA internal headers from the same OpenMM version")
    endif()
    target_sources(OpenMMColvars PRIVATE src/CudaColvarsKernel.cpp
        "${COLVARS_CORE_DIR}/cuda/colvaratoms_kernel.cu"
        "${COLVARS_CORE_DIR}/cuda/colvar_rotation_derivative_kernel.cu"
        "${COLVARS_CORE_DIR}/cuda/colvartypes_kernel.cu"
        "${COLVARS_CORE_DIR}/cuda/colvarcomp_distance_kernel.cu"
        "${COLVARS_CORE_DIR}/cuda/colvarcomp_coordnums_kernel.cu")
    target_compile_definitions(OpenMMColvars PUBLIC COLVARS_CUDA OPENMM_COLVARS_ENABLE_CUDA)
    target_include_directories(OpenMMColvars PRIVATE "${COLVARS_CORE_DIR}/cuda")
    target_link_libraries(OpenMMColvars PUBLIC CUDA::cudart "${OPENMM_CUDA_LIBRARY}")
    if(TARGET CUDA::nvtx3)
        target_link_libraries(OpenMMColvars PUBLIC CUDA::nvtx3)
    endif()
    set_property(TARGET OpenMMColvars PROPERTY CUDA_STANDARD 17)
    set_property(TARGET OpenMMColvars PROPERTY CUDA_STANDARD_REQUIRED ON)
endif()

if(OPENMM_COLVARS_BUILD_PYTHON)
    find_program(SWIG_EXECUTABLE swig REQUIRED)
    set(PYTHON_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/python")
    file(MAKE_DIRECTORY "${PYTHON_OUTPUT_DIR}")
    add_custom_command(OUTPUT "${PYTHON_OUTPUT_DIR}/openmm_colvars_wrap.cpp"
        "${PYTHON_OUTPUT_DIR}/openmm_colvars.py"
        COMMAND "${SWIG_EXECUTABLE}" -c++ -python
            "-I${CMAKE_CURRENT_SOURCE_DIR}/src" "-I${CMAKE_CURRENT_SOURCE_DIR}/python"
            -outdir "${PYTHON_OUTPUT_DIR}" -o "${PYTHON_OUTPUT_DIR}/openmm_colvars_wrap.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/python/openmm_colvars.i"
        DEPENDS python/openmm_colvars.i python/openmm_colvars_openmm_import.i src/ColvarsForce.h
        VERBATIM)
    Python_add_library(_openmm_colvars MODULE "${PYTHON_OUTPUT_DIR}/openmm_colvars_wrap.cpp")
    target_link_libraries(_openmm_colvars PRIVATE OpenMMColvars Python::Module)
    set_target_properties(_openmm_colvars PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${PYTHON_OUTPUT_DIR}"
        BUILD_RPATH "${OPENMM_LIBRARY_DIR};${OPENMM_LIBRARY_DIR}/plugins"
        INSTALL_RPATH "${OPENMM_LIBRARY_DIR};${OPENMM_LIBRARY_DIR}/plugins")
    if(NOT OPENMM_COLVARS_PYTHON_INSTALL_DIR)
        execute_process(COMMAND "${Python_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib'))"
            OUTPUT_VARIABLE OPENMM_COLVARS_PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
            COMMAND_ERROR_IS_FATAL ANY)
    endif()
    install(TARGETS _openmm_colvars LIBRARY DESTINATION "${OPENMM_COLVARS_PYTHON_INSTALL_DIR}")
    install(FILES "${PYTHON_OUTPUT_DIR}/openmm_colvars.py" DESTINATION "${OPENMM_COLVARS_PYTHON_INSTALL_DIR}")
endif()
