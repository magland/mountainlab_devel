#!/bin/bash

path0=output_r3
mountainsort sort.msh $path0/pre0.mda $path0
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
mountainsort view.msh $path0
