#!/bin/bash

COLVARS_WEBSITE_TREE=${1}
if [ -z "${COLVARS_WEBSITE_TREE}" ] ; then
    echo "Error: need the full path to the website repository as argument." >&2
    exit 1
fi

if grep -sq colvars ${COLVARS_WEBSITE_TREE}/index.html ; then
    export COLVARS_WEBSITE_TREE
else
    echo "Error: ${COLVARS_WEBSITE_TREE} does not seem to contain a Colvars webpage repository." >&2
    exit 1
fi

# If needed, initialize COLVARS_RELEASE in the same way as Makefile would
if [ -z "${COLVARS_RELEASE}" ] ; then
    COLVARS_RELEASE=$(git symbolic-ref --short -q HEAD)
fi
if [ -z "${COLVARS_RELEASE}" ] ; then
    # If we are not working on a branch, try a tag instead
    COLVARS_RELEASE=$(git describe --tags --exact-match)
fi
if [ -z "${COLVARS_RELEASE}" ] ; then
    echo "Error: ${COLVARS_RELEASE} is undefined and could not be configured from a Git branch or tag." >&2
    exit 1
fi

if make install && pushd ${COLVARS_WEBSITE_TREE} ; then
    git add ${COLVARS_RELEASE}
    if [ "x${COLVARS_RELEASE}" == "xmaster" ] ; then
        git add doxygen
        git add colvars-refman-*
    fi
    git commit -m "Update doc for version \"${COLVARS_RELEASE}\""
    # Hard-code the website repository
    git pull --rebase git@github.com:Colvars/colvars.github.io master
    git push git@github.com:Colvars/colvars.github.io master
    popd
fi
