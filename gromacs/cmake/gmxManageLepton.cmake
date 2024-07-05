# Add Lepton library, which is developed and distributed as part of OpenMM:
# https://github.com/openmm/openmm

gmx_option_multichoice(GMX_USE_LEPTON
    "Build the Lepton library interfaced with GROMACS"
    INTERNAL
    INTERNAL NONE)
mark_as_advanced(GMX_USE_LEPTON)

function(gmx_manage_lepton)
  if(GMX_USE_LEPTON STREQUAL "INTERNAL")

    set(LEPTON_DIR "${CMAKE_SOURCE_DIR}/src/external/lepton")

    file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/*.cpp)
    add_library(lepton_objlib OBJECT ${LEPTON_SOURCES})
    target_include_directories(lepton_objlib PRIVATE ${LEPTON_DIR}/include)
    # TODO support building as shared library
    target_compile_options(lepton_objlib PRIVATE -DLEPTON_BUILDING_STATIC_LIBRARY)
    set_target_properties(lepton_objlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

    add_library(lepton INTERFACE)
    target_sources(lepton INTERFACE $<TARGET_OBJECTS:lepton_objlib>)
    target_include_directories(lepton SYSTEM INTERFACE $<BUILD_INTERFACE:${LEPTON_DIR}>)

    # Set flags so that Colvars can leverage Lepton functionality
    # TODO handle the case when Lepton is built without Colvars?
    target_include_directories(colvars_objlib PRIVATE ${LEPTON_DIR}/include)
    target_compile_options(colvars_objlib PRIVATE -DLEPTON -DLEPTON_USE_STATIC_LIBRARIES)

  else()

    # Dummy target
    add_library(lepton INTERFACE)
  endif()
endfunction()
