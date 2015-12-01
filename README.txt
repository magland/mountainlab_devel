= Usage Instructions =

**** Quick Start for Ubuntu Users ******************

> sudo apt-get install git libfftw3-dev

git clone -b dev-11-14-15 https://github.com/magland/mountainsort.git
> cd mountainsort

To install Qt5, run:
./install_qt_ubuntu.sh

To compile the software run:
./compile_mountainview.sh
./compile_mountainsort_cpp.sh
./compile_spikespy.sh

Copy raw data ms11d45.dat to example_data/

Open MATLAB
> cd mountainsort
> addpath example_mscmd
> example_mscmd

(Should take less than 10 minutes, will create data in example_data)

> example_mscmd_gview

****************************************************

