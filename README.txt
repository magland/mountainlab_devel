= Usage Instructions =

**** Quick Start for Ubuntu Users ******************

git clone https://github.com/magland/mountainsort.git
> cd mountainsort
Run ./install_qt_ubuntu.sh to install Qt5
Run ./compile_mountainview_ubuntu.sh to compile mountainview
Run ./download_data_ubuntu to download some raw data
From MATLAB:
> cd mountainsort
> addpath example1
> ms_example1_process
(Should take less than 10 minutes, will create data in example1_output)
> ms_example1_view
(The mountainview program should open, see below)

****************************************************

== General Info ==

Prerequisites:
* MATLAB
* Qt4 or Qt5 (for viewing)
* Linux or Mac (has been tested on Ubuntu)

Overview:
Run the processing using ms_example1_process
View the results using ms_example1_view
You can also read the output arrays from MATLAB
The data are in example1_output

== Run mountainsort for the first time ==

From MATLAB:
> cd mountainsort
> addpath example1
> ms_example1_process

You will need to download some data into example1_data
Take a look at download_data_matlab

== Compile and run mountainview ==

You need Qt4 or Qt5 to compile mountainview.

For Ubuntu 14.04 (and perhaps other versions) use
sudo apt-add-repository ppa:ubuntu-sdk-team/ppa
sudo apt-get update
sudo apt-get install qtdeclarative5-dev
sudo apt-get install qt5base-dev qtscript5-dev make g++

Otherwise, download/install Qt5 from here:
http://www.qt.io/download/
You should use the open source version

Once installed, qmake should be in your path.
Test that by running "qmake -version"
If not, then contact Jeremy for help.

From console:
> cd mountainsort/mountainview/src
> qmake
> make

This should create an executable:
mountainsort/mountainview/bin/mountainview

From MATLAB:
> cd mountainsort
> addpath example1
> example1_mountainview

You should get a small window with a few buttons.
Press the buttons to view different aspects of the data

== Accessing the output data ==

The results will be in ms11d45A/output
Use util/readmda.m and util/writemda.m to
read and write the .mda file format.
See http://magland.github.io//articles/mda-format/

raw.mda: the raw, pre-processed data
times.mda/labels.mda: the detected event times and corresponding neuron labels
locations.mda: geometric locations of the electrodes
adjacency.mda: adjacency matrix for the electrodes - obtained from locations.mda
templates.mda: the templates for each spike type
primary_channels.mda: the primary channels for each spike template
spikespy.raw.* folder: a temporary folder created by the spikespy viewer
cross_correlograms.mda: data for the cross-correlograms -- see example1_view_cross_correlograms.m







