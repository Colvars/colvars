if(BUILD_TESTS)

  # Build functional tests' executables
  add_subdirectory(${COLVARS_SOURCE_DIR}/tests/functional tests/functional)

  if(DEFINED CMAKE_SYSTEM_NAME)

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")

      # Create a separate target for copying test files
      add_custom_target(copy_test_files ALL
        COMMENT "Copying test input files and configurations"
        VERBATIM
      )

      # TODO automate this step and make it OS-portable
      # add_custom_command(TARGET colvars POST_BUILD
      #   COMMAND bash build_tests.sh
      #   ${CMAKE_BINARY_DIR}/tests/functional
      #   WORKING_DIRECTORY ${COLVARS_SOURCE_DIR}/tests)

      # Copy the Colvars configuration files
      file(GLOB TEST_CONFIG_FILES ${COLVARS_SOURCE_DIR}/tests/input_files/*/test.in)
      if(NOT COLVARS_TORCH)
        # TODO create a way to detect test dependencies at some point
        list(REMOVE_ITEM TEST_CONFIG_FILES ${COLVARS_SOURCE_DIR}/tests/input_files/torchann-dihedral_harmonic-fixed/test.in)
      endif()

      foreach(TEST_CONFIG_FILE ${TEST_CONFIG_FILES})
        get_filename_component(TEST_NAME ${TEST_CONFIG_FILE} DIRECTORY)
        get_filename_component(TEST_NAME ${TEST_NAME} NAME)
        add_test(NAME ${TEST_NAME}
          COMMAND run_colvars_test ${TEST_NAME}/test.in trajectory.xyz "${TEST_NAME}/test_out"
          WORKING_DIRECTORY
          ${CMAKE_CURRENT_BINARY_DIR}/tests/functional)
        add_test(NAME "${TEST_NAME}_spiff"
          COMMAND bash compare_test.sh ${TEST_NAME}
          WORKING_DIRECTORY
          ${CMAKE_CURRENT_BINARY_DIR}/tests/functional)
      endforeach()

      # Copy input files (coordinates, index files, etc)
      add_custom_command(TARGET copy_test_files
        COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${COLVARS_SOURCE_DIR}/tests/input_files
          ${CMAKE_BINARY_DIR}/tests/functional
        COMMENT "Copying test input files"
      )

      # Make the main target depend on the copy target
      add_dependencies(colvars copy_test_files)

    endif()
  endif()
  enable_testing()
endif()

if(BUILD_UNITTESTS)
  # Build unit tests executables
  add_subdirectory(${COLVARS_SOURCE_DIR}/tests/unittests tests/unittests)
endif()
