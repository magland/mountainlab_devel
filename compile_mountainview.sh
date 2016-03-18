#!/bin/bash

# Compile mountainview
echo "Compiling mountainview"
cd mountainlab/mountainview/src
qmake
make -j 8
cd ../..


