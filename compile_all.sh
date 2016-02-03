#!/bin/bash

MY_PATH="`dirname \"$BASH_SOURCE\"`"
$MY_PATH/compile_mountainsort_cpp.sh
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
$MY_PATH/compile_mountainview.sh
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
$MY_PATH/compile_spikespy.sh
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi


