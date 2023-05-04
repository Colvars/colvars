# General build script suitable for cross-platform CI tests

if(DEFINED CMAKE_SCRIPT_MODE_FILE)
  get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_SCRIPT_MODE_FILE} DIRECTORY)
  get_filename_component(COLVARS_SOURCE_DIR ${COLVARS_SOURCE_DIR} DIRECTORY)
else()
  get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
endif()

set(COLVARS_LEPTON ON)
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
else()
  if(${CMAKE_CXX_STANDARD} GREATER 70)
    set(COLVARS_LEPTON OFF)
  endif()
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
  # We need a shared library to load the Tcl package
  set(BUILD_SHARED_LIBS ON)
endif()

if(NOT DEFINED COLVARS_TCL)
  set(COLVARS_TCL ON)
endif()

if(DEFINED CMAKE_SYSTEM_NAME AND COLVARS_TCL)

  # If available, use pre-downloaded OS-specific TCL libraries

  if(EXISTS "${COLVARS_SOURCE_DIR}/devel-tools/packages")
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
      set(TCL_DIR "${COLVARS_SOURCE_DIR}/devel-tools/packages/tcl8.5.9-linux-x86_64")
      set(TCL_LIBRARY "libtcl8.5.a")
    endif()

    if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
      set(TCL_DIR "${COLVARS_SOURCE_DIR}/devel-tools/packages/tcl8.5.9-win64")
      set(TCL_LIBRARY "tcl85.lib")
    endif()

    if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
      set(TCL_DIR "${COLVARS_SOURCE_DIR}/devel-tools/packages/tcl8.5.9-macosx-x86_64-threaded")
      set(TCL_LIBRARY "libtcl8.5.a")
    endif()

    set(DEFINE_TCL_DIR "-DTCL_DIR=${TCL_DIR}")
    set(DEFINE_TCL_LIBRARY "-DTCL_LIBRARY=${TCL_DIR}/lib/${TCL_LIBRARY}")
  endif()

endif()

if(COLVARS_LEPTON)
  if(NOT DEFINED LEPTON_DIR)
    set(LEPTON_DIR "${COLVARS_SOURCE_DIR}/openmm-source/libraries/lepton")
    if(NOT EXISTS ${LEPTON_DIR})
      execute_process(COMMAND git clone --depth=1 https://github.com/openmm/openmm.git "${COLVARS_SOURCE_DIR}/openmm-source")
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
  set(DEFINE_CC_CCACHE "-DCMAKE_C_COMPILER_LAUNCHER=$ENV{CCACHE}")
  set(DEFINE_CXX_CCACHE "-DCMAKE_CXX_COMPILER_LAUNCHER=$ENV{CCACHE}")
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND}
  -S cmake
  -B ${BUILD_DIR}
  -D COLVARS_DEBUG=${COLVARS_DEBUG}
  -D CMAKE_BUILD_TYPE=RelWithDebinfo
  -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
  -D WARNINGS_ARE_ERRORS=ON
  -D CMAKE_VERBOSE_MAKEFILE=ON
  -D ENABLE_COVERAGE=ON
  -D CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
  ${DEFINE_CC_CCACHE}
  ${DEFINE_CXX_CCACHE}
  -D COLVARS_TCL=${COLVARS_TCL}
  ${DEFINE_TCL_DIR}
  ${DEFINE_TCL_LIBRARY}
  -D COLVARS_LEPTON=${COLVARS_LEPTON}
  -D LEPTON_DIR=${LEPTON_DIR}
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

if(NOT result EQUAL 0)
  message(FATAL_ERROR "Error building library.")
else()
  execute_process(
    COMMAND ${CMAKE_CTEST_COMMAND}
    --test-dir ${BUILD_DIR}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Error running library tests.")
  endif()
endif()
