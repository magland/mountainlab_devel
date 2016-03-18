#!/bin/bash

path0=output_ms11d45
mountainsort ds001_sort.msh $path0/pre0.mda $path0 --adjacency_matrix=$path0/adjacency_matrix.mda
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
mountainsort ds001_view.msh $path0
