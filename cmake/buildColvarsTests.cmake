if(BUILD_TESTS)
  add_subdirectory(${COLVARS_SOURCE_DIR}/tests/functional tests/functional)
  if(DEFINED CMAKE_SYSTEM_NAME)
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
      message(STATUS "Copying input files for tests")
      # TODO make this portable outside Linux
      add_custom_command(TARGET colvars POST_BUILD
        COMMAND bash build_tests.sh
        ${CMAKE_BINARY_DIR}/tests/functional
        WORKING_DIRECTORY ${COLVARS_SOURCE_DIR}/tests)
      file(GLOB TEST_INPUT_FILES ${COLVARS_SOURCE_DIR}/tests/input_files/*)
      message("${TEST_INPUT_FILES}")
      foreach(TEST_INPUT_FILE ${TEST_INPUT_FILES})
        add_custom_command(TARGET colvars POST_BUILD
          COMMAND ${CMAKE_COMMAND}
          -E copy ${TEST_INPUT_FILE} ${CMAKE_BINARY_DIR}/tests/functional)
      endforeach(TEST_INPUT_FILE)
    endif()
  endif()
endif()

option(BUILD_UNITTESTS "Build unit tests" ${BUILD_TESTS})
if(BUILD_UNITTESTS)
  add_subdirectory(${COLVARS_SOURCE_DIR}/tests/unittests tests/unittests)
endif()
