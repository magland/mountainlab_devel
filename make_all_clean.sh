#!/bin/bash

echo "Cleaning all three Qt projects back to source files"

(cd mountainview/src/spikespy/src; make clean)
(cd mountainview/src; make clean)
(cd cpp/src; make clean)
