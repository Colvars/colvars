if(COLVARS_LEPTON)

  if(NOT EXISTS ${LEPTON_DIR})
    message(FATAL_ERROR "With -DCOLVARS_LEPTON=ON, the lepton folder must be copied into ${COLVARS_SOURCE_DIR} or provided by -DLEPTON_DIR=xxx")
  endif()

  if(${CMAKE_CXX_STANDARD} GREATER 70) # Yet another Y2K idiosyncrasy
    message(FATAL_ERROR "With -DCOLVARS_LEPTON=ON, CMAKE_CXX_STANDARD must be 11 or later")
  endif()

  file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
  add_library(lepton OBJECT ${LEPTON_SOURCES})
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(lepton PRIVATE /DLEPTON_BUILDING_STATIC_LIBRARY /wd4267)
  endif()
  set_property(TARGET lepton PROPERTY POSITION_INDEPENDENT_CODE 1)
  target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(colvars PRIVATE /DLEPTON /DLEPTON_USE_STATIC_LIBRARIES)
  else()
    target_compile_options(colvars PRIVATE -DLEPTON)
  endif()
  target_include_directories(colvars PRIVATE ${LEPTON_DIR}/include)
  target_link_libraries(colvars lepton)

  # Silence warnings for Lepton
  target_compile_options(lepton PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Wno-tautological-undefined-compare -Wno-unknown-warning-option>)

endif()