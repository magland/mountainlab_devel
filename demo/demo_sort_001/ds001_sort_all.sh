#!/bin/bash

path0=output_tetrode1
mountainsort ds001_sort.msh $path0/pre0.mda $path0 --detectability_threshold=4
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
mountainsort ds001_view.msh $path0

path0=output_tetrode2
mountainsort ds001_sort.msh $path0/pre0.mda $path0 --detectability_threshold=4
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
mountainsort ds001_view.msh $path0

path0=output_ms11d45
mountainsort ds001_sort_ms11d45.msh $path0/pre0.mda $path0 --detectability_threshold=4
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
mountainsort ds001_view.msh $path0

path0=output_synth_EJ_K7
mountainsort ds001_sort.msh $path0/pre0.mda $path0 --detectability_threshold=5
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
mountainsort ds001_view.msh $path0
