if(BUILD_TESTS)

  # Build functional tests' executables
  add_subdirectory(${COLVARS_SOURCE_DIR}/tests/functional tests/functional)

  if(DEFINED CMAKE_SYSTEM_NAME)

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")

      # TODO automate this step and make it OS-portable
      # add_custom_command(TARGET colvars POST_BUILD
      #   COMMAND bash build_tests.sh
      #   ${CMAKE_BINARY_DIR}/tests/functional
      #   WORKING_DIRECTORY ${COLVARS_SOURCE_DIR}/tests)

      # Copy the Colvars configuration files
      file(GLOB TEST_CONFIG_FILES ${COLVARS_SOURCE_DIR}/tests/input_files/*/test.in)
      foreach(TEST_CONFIG_FILE ${TEST_CONFIG_FILES})
        get_filename_component(TEST_NAME ${TEST_CONFIG_FILE} DIRECTORY)
        get_filename_component(TEST_NAME ${TEST_NAME} NAME)
        add_custom_command(TARGET colvars POST_BUILD
          COMMAND ${CMAKE_COMMAND}
          -E copy ${TEST_CONFIG_FILE}
          ${CMAKE_BINARY_DIR}/tests/functional/${TEST_NAME}/test.in)
        add_test(NAME ${TEST_NAME}
          COMMAND run_colvars_test ${TEST_NAME}/test.in
          WORKING_DIRECTORY
          ${CMAKE_CURRENT_BINARY_DIR}/tests/functional)
      endforeach()

      # Copy other input files (coordinates, index files, etc)
      file(GLOB TEST_INPUT_FILES ${COLVARS_SOURCE_DIR}/tests/input_files/*)
      foreach(TEST_INPUT_FILE ${TEST_INPUT_FILES})
        add_custom_command(TARGET colvars POST_BUILD
          COMMAND ${CMAKE_COMMAND}
          -E copy ${TEST_INPUT_FILE}
          ${CMAKE_BINARY_DIR}/tests/functional)
      endforeach(TEST_INPUT_FILE)

    endif()
  endif()
  enable_testing()
endif()

if(BUILD_UNITTESTS)
  # Build unit tests executables
  add_subdirectory(${COLVARS_SOURCE_DIR}/tests/unittests tests/unittests)
endif()
