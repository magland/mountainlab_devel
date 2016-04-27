# mountainlab_devel: spike-sorting validation and development codes for mountainlab

Authors: Jeremy Magland and Alex Barnett

Started February 2016. Version 4/11/16

See [`mountainlab`](https://github.com/magland/mountainlab)  

This is a collection of MATLAB codes for running tests on `mountainlab`
(abbreviated by ML from now on),
which is a C++/qt5-based package for automatic spike-sorting of extracellular
neuronal recordings, and visualization.

Currently accuracy is measured on various synthetic and intracellular
ground-truthed datasets.
Stability-based validation is not yet implemented; for now see 
[validspike](https://github.com/ahbarnett/validspike) for this.

### Requirements

* recent MATLAB (eg R2012b or newer). No toolboxes are needed.  

* `mountainlab`. Be sure to compile mountainsort and mountainview.


### Installation

As usual, `git clone https://github.com/magland/mountainlab_devel`

Make sure you have `mountainlab` installed and compiled. Optionally
install `spikespy` in at `~/spikespy`

Start MATLAB from the `mountainlab_devel` directory.

Make sure `mountainlab/matlab` is in your path, and run
`mountainlab_setup`

Run `ms_setup`

To check the MATLAB code works, run `simplesorter`. If spikespy
is installed you'll see a window pop up.

To check the interface to ML works, run eg `accuracy_jfm_april_sort`


### Guide to directories

`unit_tests` : basic tests of various ML components, and generation of the demo (default) timeseries dataset. Warning: this creates the directories at and below `unit_tests/demo_data` and fills them with a few tens of megabytes of files. (Also many of the `ms_*` commands in ML do self-tests when run without arguments).  

`unit_tests/testall.m` is a useful list of all MATLAB-driven self tests.  

`examples` : examples of drivers for accuracy tests and spike sorters.  

`sorting_algs` : MATLAB wrappers for various sorting algorithms in a standard interface.  

`validation` : functions needed for accuracy tests and comparing sorter outputs.   

`view` : various MATLAB-based viewer functions independent of ML.  


Test codes also write to the location given by MATLAB's `tempdir` command,
which is usually `\tmp` in a linux system.


### Setting up external datasets

A synthetic data generator is shipped with waveforms deriving from
a tiny retinal dataset from 2005 due to E J Chichilnisky's group.
However, validation can
also be done against various external datasets which are not supplied
in this repo. Here's how to set them up:

Create an external dataset directory which we'll call extdata. From
the `mountainlab_devel` directory do

`ln -s extdata ext_datasets`  

If you are on the SCDA network you may simply replace
extdata above with `~ahb/ss_datasets` and be done.

If you are anywhere else and want each of the datasets, do the following:


#### Harris 2000 tetrode, 4 minutes, 4 EC channels, 1 IC channel (38 MB)

Get a [CRCNS account](https://crcns.org/)

Download from [here](https://crcns.org/data-sets/hc/hc-1)
the file `d5331.zip` to
`extdata/Harris2000/d5331/`
Extract it. This creates amongst other files,
`extdata/Harris2000/d5331/d533101.dat`
which is the only one used.

Info is [here](https://crcns.org/data-sets/hc/hc-1/about)

Test with `simplesorter_harris2000`


#### Martinez et al 2009, 2 minutes, 1 EC channel, simulated (5 x 22 MB)

Create `extdata/MartinezQuiroga2009_sims` and into download via

`wget http://www2.le.ac.uk/departments/engineering/research/bioengineering/neuroengineering-lab/simulations/simulation-n.mat`

where `n` is 1,..,5.

Test with `simplesorter_martinez2009`.  The argument to `grab_martinez2009_dataset` is `n`.


#### NeuroCube (Camunas-Mesa et al 2013), tetrode, simulated (12 MB, 115 MB)

Creation of datasets:
Download [`NeuroCube'](http://www.le.ac.uk/csn/neurocube)
MATLAB package from [Leicester Univ](http://www.le.ac.uk/large_dl/NeuroCube.zip).
(Luis Camunas-Mesa also provided versions of `exprnd` and `gprnd` that remove the need for the
Stats Toolbox. The Parallel Toolbox is made good use of if available.)
Unpack, go into `mat` directory, run `neurocube`, choose Tetrode and leave all else default.
This chooses 300000 neurons/mm^3, 7% rate of active neurons, exponential pdf of firing rates.
Click Generate cube, then Run (takes a few minutes), then Save as `Spikesim_default_tet_4-25-16.mat`;
this is a 30 second dataset with around 600 ground-truth neurons (total 20000 neurons simulated).
To make a 5-minute dataset, enter 300 in Duration, and rerun (run-time around 40 minutes, even with Parallel Toolbox), Save as `Spikesim_tet_5min.mat`.
Place these files in `extdata/NeuroCube/mat/`; this has been done on the SCDA network.

`grab_martinez2009_dataset(n)` with `n=1` accesses the 30-second dataset, and `n=2` the 5-minute dataset.

Test with `accuracyneurocube_allsorters`



### To do list

* more ground-truthed datasets

* stability metrics

