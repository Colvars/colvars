#!/bin/bash

# Get the top Colvars repo directory
TOPDIR=$(git rev-parse --show-toplevel)
if [ ! -d ${TOPDIR} ] ; then
  echo "Error: cannot identify top project directory." >& 2
  exit 1
fi

# Get the build directory
COLVARS_BUILD_DIR="$(pwd)"

# Get the clang-tidy binary
if [[ -v CLANG_TIDY ]]; then
  CLANG_TIDY_BINARY=$CLANG_TIDY
  echo "Using $CLANG_TIDY_BINARY as clang-tidy"
else
  CLANG_TIDY_BINARY="clang-tidy"
fi

# Get the Colvars source directory
COLVARS_SOURCE_DIR="$TOPDIR/src"
# Find all source files that are required to build the run_colvars_test
COLVARS_STUB_SOURCE_DIR="$TOPDIR/misc_interfaces/stubs"
COLVARS_TEST_SOURCE_DIR="$TOPDIR/tests/functional/"
COLVARS_SOURCE_FILES=($(find $COLVARS_STUB_SOURCE_DIR $COLVARS_TEST_SOURCE_DIR $COLVARS_SOURCE_DIR -name "*.cpp"))

printf '%s\n' ${COLVARS_SOURCE_FILES[@]} | \
    xargs -I{} -P $(nproc) \
    bash -c "$CLANG_TIDY_BINARY -p=$COLVARS_BUILD_DIR -header-filter=.* --warnings-as-errors='*' '{}' > '{}'.log"

retcode=$?

if [ $retcode != 0 ] ; then
    cat "${COLVARS_SOURCE_FILES[@]/%/.log}"
    exit $retcode
fi

exit 0
