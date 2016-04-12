# mountainlab_devel: spike-sorting validation and development codes for mountainlab

Authors: Jeremy Magland and Alex Barnett

Started February 2016. Version 4/11/16

See [`mountainlab`](https://github.com/magland/mountainlab)  

This is a collection of MATLAB codes for running tests on `mountainlab`
(abbreviated by ML from now on),
which is a C++/qt5-based package for automatic spike-sorting of extracellular
neuronal recordings, and visualization.

Currently accuracy is measured.
Stability-based validation is not yet implemented; see for this
[validspike](https://github.com/ahbarnett/validspike)

### Requirements

* recent MATLAB (eg R2012b or newer). No toolboxes are needed.  

* `mountainlab`. Be sure to compile mountainsort and mountainview.


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


### To do list

* more ground-truthed datasets

* stability metrics

