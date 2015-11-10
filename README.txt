= Usage Instructions =

Prerequisites:
* MATLAB
* Qt4 or Qt5 (for viewing)
* Linux or Mac (has been tested on Ubuntu)

Overview:
Run the processing using mountainsort_example
View the results using mountainview_example
Read the output arrays from MATLAB


== Run mountainsort for the first time ==

From MATLAB:
> cd mountainsort
> mountainview_example1

You will get a message telling you to download some data.
Follow those instructions and then try re-running.

Processing should take less than one minute


== Compile and run mountainview ==

You need Qt4 or Qt5 to compile mountainview.

Download/install Qt5 from here:
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
> mountainview_example1

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







