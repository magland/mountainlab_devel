#!/bin/bash

path0=output_tetrode2
mountainsort ds002_sort.msh $path0/pre0.mda $path0
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
mountainsort ds002_view.msh $path0
