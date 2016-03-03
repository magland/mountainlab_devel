#!/bin/bash

path0=output_synth_EJ_K7
mountainsort ds001_sort.msh $path0/pre0.mda $path0 --detectability_threshold=5
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
mountainsort ds001_view.msh $path0
