add_executable(colvarvalue_unit3vector colvarvalue_unit3vector.cpp)
target_link_libraries(colvarvalue_unit3vector PRIVATE colvars)
target_include_directories(colvarvalue_unit3vector PRIVATE ${COLVARS_SOURCE_DIR}/src)
add_test(NAME colvarvalue_unit3vector COMMAND colvarvalue_unit3vector)

add_executable(file_io file_io.cpp)
target_link_libraries(file_io PRIVATE colvars)
target_include_directories(file_io PRIVATE ${COLVARS_SOURCE_DIR}/src)
add_test(NAME file_io COMMAND file_io)

if(COLVARS_TCL)
  add_executable(embedded_tcl embedded_tcl.cpp)
  target_link_libraries(embedded_tcl PRIVATE colvars)
  target_include_directories(embedded_tcl PRIVATE ${COLVARS_SOURCE_DIR}/src)
  add_test(NAME embedded_tcl COMMAND embedded_tcl)
endif()

# Last test requires the stubs proxy
add_executable(read_xyz_traj read_xyz_traj.cpp)
set(COLVARS_STUBS_DIR ${COLVARS_SOURCE_DIR}/tests/stubs/)
target_link_libraries(read_xyz_traj PRIVATE colvars colvars_stubs)
target_include_directories(read_xyz_traj PRIVATE ${COLVARS_SOURCE_DIR}/src)
target_include_directories(read_xyz_traj PRIVATE ${COLVARS_STUBS_DIR})
add_test(NAME read_xyz_traj COMMAND read_xyz_traj)
add_custom_command(
        TARGET read_xyz_traj POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/da-traj.xyz
                ${CMAKE_CURRENT_BINARY_DIR}/da-traj.xyz)