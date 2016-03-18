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

Copy raw data to franklab/raw

Open MATLAB
> cd mountainsort
> ms_setup_path
> addpath franklab/[name_of_experiment]
> [name_of_experiment]

(Should take less than 3 minutes, will create data in franklab/[name_of_experiment]/output*)

****************************************************

= Examples =

examples  directory contains these

= Unit tests =

All commands ms_*.m do self-test when called without arguments
test_filter: compares MATLAB and C filter pipelines

= Validation = 

ln -s your_external_dataset_dir ext_datasets
