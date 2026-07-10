#!/bin/bash
# This script generates reference output for the AutoDiff tests.
# Run from the tests/input_files directory, with the test names as arguments, e.g.:
#   ./generate_ref.sh test1 test2 test3

RUN=~/Projects/colvars/build/tests/functional/run_colvars_test
TESTS=$*

for TEST in $TESTS
do
    # remove trailing _spiff if present (output from ctest)
    TEST=${TEST%_spiff}
    $RUN -c $TEST/test.in  -t trajectory.xyz --force -o $TEST/AutoDiff/test_out
done
