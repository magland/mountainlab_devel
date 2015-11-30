#!/bin/bash

# Compile spikespy
echo "Compiling spikespy"
cd mountainview/src/spikespy/src
qmake
make -j 8
cd ../../../..


