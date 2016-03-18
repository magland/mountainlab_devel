#!/bin/bash

# Compile spikespy
echo "Compiling spikespy"
cd OBSOLETE_mountainview/src/spikespy/src
qmake
make -j 8
cd ../../../..


