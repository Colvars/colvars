#!/bin/bash

if [ ! -d Common ] ; then
  ln -s ../library/Common
fi

../library/run_tests.sh "$@"


