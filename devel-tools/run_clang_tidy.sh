#!/bin/bash

# Get the top Colvars repo directory
TOPDIR=$(git rev-parse --show-toplevel)
if [ ! -d ${TOPDIR} ] ; then
  echo "Error: cannot identify top project directory." >& 2
  exit 1
fi

# Get the build directory
COLVARS_BUILD_DIR="$(pwd)"

# Get the Colvars source directory
COLVARS_SOURCE_DIR="$TOPDIR/src"
# Find all source files that are required to build the run_colvars_test
COLVARS_STUB_SOURCE_DIR="$TOPDIR/misc_interfaces/stubs"
COLVARS_TEST_SOURCE_DIR="$TOPDIR/tests/functional/"
COLVARS_SOURCE_FILES=$(find $COLVARS_STUB_SOURCE_DIR $COLVARS_TEST_SOURCE_DIR $COLVARS_SOURCE_DIR -name "*.cpp")

num_test_failed=0
all_output=""
# Run clang-tidy over all files and save the results
for colvars_src_file in $COLVARS_SOURCE_FILES; do
  clang_tidy_command="clang-tidy -p=$COLVARS_BUILD_DIR --warnings-as-errors='*' $colvars_src_file"
  echo "Running $clang_tidy_command"
  output="$(eval $clang_tidy_command 2>&1)" || exit_code=$?
  all_output+="$output"
  if [[ $exit_code -gt 0 ]]; then
    num_test_failed=`expr $num_test_failed + 1`
  fi
done

if [[ $num_test_failed -gt 0 ]]; then
  echo "There are $num_test_failed test(s) failed."
  echo "$all_output"
  exit 1
else
  exit 0
fi
