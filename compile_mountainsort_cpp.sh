#!/bin/bash

# Compile mountainsort/cpp
echo "Compiling mountainsort/cpp"
cd cpp/src
qmake
make -j 8
cd ../..


