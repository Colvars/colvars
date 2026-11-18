set(COLVARS_STUBS_DIR ${COLVARS_SOURCE_DIR}/misc_interfaces/stubs/)

foreach(CMD
    colvarvalue_unit3vector
    file_io
    memory_stream
    read_xyz_traj
    parse_error
    triclinic_minimum_image
  )
  add_executable(${CMD} ${CMD}.cpp)
  target_link_libraries(${CMD} PRIVATE colvars)
  target_include_directories(${CMD} PRIVATE ${COLVARS_SOURCE_DIR}/src)
  target_include_directories(${CMD} PRIVATE ${COLVARS_STUBS_DIR})
  target_link_libraries(${CMD} PRIVATE colvars colvars_stubs)
  add_test(NAME ${CMD} COMMAND ${CMD})
endforeach()

if(COLVARS_TCL)
  add_executable(embedded_tcl embedded_tcl.cpp)
  target_link_libraries(embedded_tcl PRIVATE colvars)
  target_include_directories(embedded_tcl PRIVATE ${COLVARS_SOURCE_DIR}/src)
  add_test(NAME embedded_tcl COMMAND embedded_tcl)
endif()

foreach(input IN ITEMS index.ndx test.in da-traj.xyz test_error.in)
  configure_file(${input} ${input} COPYONLY)
endforeach()
