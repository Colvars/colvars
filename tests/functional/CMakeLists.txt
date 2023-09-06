set(COLVARS_STUBS_DIR ${COLVARS_SOURCE_DIR}/misc_interfaces/stubs/)
add_library(colvars_stubs OBJECT ${COLVARS_STUBS_DIR}/colvarproxy_stub.cpp)
target_include_directories(colvars_stubs PRIVATE ${COLVARS_SOURCE_DIR}/src)


add_executable(run_colvars_test run_colvars_test.cpp)
target_link_libraries(run_colvars_test PRIVATE colvars colvars_stubs)
target_include_directories(run_colvars_test PRIVATE ${COLVARS_SOURCE_DIR}/src)
target_include_directories(run_colvars_test PRIVATE ${COLVARS_STUBS_DIR})
