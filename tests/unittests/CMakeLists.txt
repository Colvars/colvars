set(COLVARS_STUBS_DIR ${COLVARS_SOURCE_DIR}/misc_interfaces/stubs/)

foreach(CMD
    colvarvalue_unit3vector
    file_io
    memory_stream
    read_xyz_traj
    parse_error
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

add_custom_command(
        TARGET read_xyz_traj POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/index.ndx
                ${CMAKE_CURRENT_BINARY_DIR}/index.ndx
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/test.in
                ${CMAKE_CURRENT_BINARY_DIR}/test.in
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/da-traj.xyz
                ${CMAKE_CURRENT_BINARY_DIR}/da-traj.xyz)

add_custom_command(
        TARGET parse_error POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/test_error.in
                ${CMAKE_CURRENT_BINARY_DIR}/test_error.in)
