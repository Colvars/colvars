# General build script suitable for cross-platform CI tests

if(DEFINED CMAKE_SCRIPT_MODE_FILE)
  get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_SCRIPT_MODE_FILE} DIRECTORY)
  get_filename_component(COLVARS_SOURCE_DIR ${COLVARS_SOURCE_DIR} DIRECTORY)
else()
  get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
endif()

if(EXISTS "/opt/libtorch")
  if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
  endif()
  if(${CMAKE_CXX_STANDARD} GREATER_EQUAL 17)
    message("Enabling Torch interface, using library at /opt/libtorch")
    set(DEFINE_TORCH "-DCOLVARS_TORCH=ON")
    set(DEFINE_TORCH_PREFIX "-DLIBTORCH_PREFIX=/opt/libtorch")
  endif()
endif()

if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebinfo)
endif()

if(TRAP_FPE)
  message("Trapping floating-point exceptions in functional tests")
  set(DEFINE_TRAP_FPE "-DCOLVARS_TRAP_FPE=ON")
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
  # We need a shared library to load the Tcl package
  set(BUILD_SHARED_LIBS ON)
endif()

if(NOT DEFINED COLVARS_TCL)
  set(COLVARS_TCL ON)
endif()

# If available, use pre-downloaded TCL libraries
if(EXISTS "${COLVARS_SOURCE_DIR}/devel-tools/packages")
  if(DEFINED CMAKE_HOST_SYSTEM_NAME)
    if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin")
      execute_process(
        COMMAND uname -m
        OUTPUT_VARIABLE CMAKE_HOST_SYSTEM_PROCESSOR
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      message("CMAKE_HOST_SYSTEM_PROCESSOR = ${CMAKE_HOST_SYSTEM_PROCESSOR}")
      if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux" AND "${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
        set(TCL_DIR "${COLVARS_SOURCE_DIR}/devel-tools/packages/tcl8.6.17-linux-x86_64-threaded")
        set(TCL_LIBRARY "libtcl8.6.a")
      endif()
      if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux" AND "${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
        set(TCL_DIR "${COLVARS_SOURCE_DIR}/devel-tools/packages/tcl8.6.17-linux-arm64-threaded")
        set(TCL_LIBRARY "libtcl8.6.a")
      endif()
      if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin" AND "${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
        set(TCL_DIR "${COLVARS_SOURCE_DIR}/devel-tools/packages/tcl8.5.9-macosx-x86_64-threaded")
        set(TCL_LIBRARY "libtcl8.5.a")
      endif()
    endif()
    # Assume Intel when using Windows
    if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Windows")
      set(TCL_DIR "${COLVARS_SOURCE_DIR}/devel-tools/packages/tcl8.5.9-win64")
      set(TCL_LIBRARY "tcl85.lib")
    endif()
  endif()
endif()

if(COLVARS_TCL AND DEFINED TCL_DIR)
  set(DEFINE_TCL_DIR "-DTCL_DIR=${TCL_DIR}")
  set(DEFINE_TCL_LIBRARY "-DTCL_LIBRARY=${TCL_DIR}/lib/${TCL_LIBRARY}")
endif()


set(COLVARS_LEPTON ON)

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
endif()

if(DEFINED ENV{CMAKE_BUILD_DIR})
  set(BUILD_DIR $ENV{CMAKE_BUILD_DIR})
else()
  if(DEFINED ENV{CXX_VERSION})
    set(BUILD_DIR build_$ENV{CXX}$ENV{CXX_VERSION}_C++${CMAKE_CXX_STANDARD})
  else()
    set(BUILD_DIR build)
  endif()
endif()

if(DEFINED ENV{CCACHE})
  set(CCACHE $ENV{CCACHE})
else()
  find_program(CCACHE "ccache")
  if(CCACHE)
    set(CCACHE "ccache")
  endif()
endif()

if(NOT ${CCACHE} STREQUAL "CCACHE-NOTFOUND")
  message(STATUS "Using ${CCACHE} as the compiler launcher")
  set(DEFINE_CC_CCACHE "-DCMAKE_C_COMPILER_LAUNCHER=${CCACHE}")
  set(DEFINE_CXX_CCACHE "-DCMAKE_CXX_COMPILER_LAUNCHER=${CCACHE}")
endif()

if(NOT DEFINED ENV{CXX})
  # Use Clang if available to catch more errors
  find_program(CLANG "clang++")
  if(CLANG)
    message(STATUS "Clang is available, using it as the default compiler")
    set(DEFINE_CXX_COMPILER "-DCMAKE_CXX_COMPILER=clang++")
  endif()
endif()

# See https://stackoverflow.com/questions/26836361/check-if-generating-a-visual-studio-solution-or-makefile-from-cmake
if(CMAKE_GENERATOR MATCHES "Visual Studio")
  # Workaround for https://gitlab.kitware.com/cmake/cmake/-/issues/27116
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /GR")
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND}
  -S cmake
  -B ${BUILD_DIR}
  -D COLVARS_DEBUG=${COLVARS_DEBUG}
  -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
  -D WARNINGS_ARE_ERRORS=ON
  -D CMAKE_VERBOSE_MAKEFILE=ON
  -D CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
  -D CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
  ${DEFINE_CXX_COMPILER}
  ${DEFINE_CC_CCACHE}
  ${DEFINE_CXX_CCACHE}
  ${DEFINE_PYTHON}
  -D COLVARS_TCL=${COLVARS_TCL}
  ${DEFINE_TCL_DIR}
  ${DEFINE_TCL_LIBRARY}
  ${DEFINE_TORCH}
  ${DEFINE_TORCH_PREFIX}
  ${DEFINE_TRAP_FPE}
  -D COLVARS_LEPTON=${COLVARS_LEPTON}
  -D LEPTON_DIR=${LEPTON_DIR}
  -D CMAKE_PREFIX_PATH="/opt/libtorch/share/cmake"
  RESULT_VARIABLE result
  )

if(NOT result EQUAL 0)
  message(FATAL_ERROR "Error generating CMake buildsystem.")
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND}
  --build ${BUILD_DIR}
  --parallel
  RESULT_VARIABLE result
  )

option(RUN_TESTS "Run library tests" ON)

if(NOT result EQUAL 0)
  message(FATAL_ERROR "Error building library.")
else()
  if(RUN_TESTS)
    execute_process(
      COMMAND ${CMAKE_CTEST_COMMAND}
      WORKING_DIRECTORY ${BUILD_DIR}
      RESULT_VARIABLE result
    )
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "Error running library tests.")
    endif()
  endif()
endif()
