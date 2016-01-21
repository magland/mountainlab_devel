= Usage Instructions =

**** Quick Start for Ubuntu Users ******************

> sudo apt-get install git libfftw3-dev

git clone https://github.com/magland/mountainsort.git
> cd mountainsort

To install Qt5, run:
./install_qt_ubuntu.sh

To compile the software run:
./compile_mountainview.sh
./compile_mountainsort_cpp.sh
./compile_spikespy.sh

Copy raw data ms11d45.dat to franklab/raw

Open MATLAB
> cd mountainsort
> addpath franklab/experiment1
> experiment1_sort

(Should take less than 3 minutes, will create data in franklab/experiment1/output)

> experiment1_view

****************************************************

