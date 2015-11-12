#!/bin/bash

basedir=$PWD

# Download data if needed
echo "Downloading raw data"
path=$basedir/example1_data/ms11d45A_pre.mda
url="http://voms.simonsfoundation.org:50013/rXcGof5b0pDR4zkEDCevBLfo95RB7/frank_ucsf_example1/ms11d45A_pre.mda"
echo curl $url > $path
curl $url > $path


