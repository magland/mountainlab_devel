Demo scripts for sorting and viewing output
Author: Jeremy Magland
Date: 1 Mar 2016

=== Step 1: Prepare raw data ===

You need to place the raw data in the appropriate locations within franklab/raw

Then run ds001_prepare_raw_data.m in MATLAB. This will create some output directories and place a pre0.mda file in each.


=== Step 2: Run the sorting ===

You can either run sorting using a MATLAB script, or using the command line (mountainsort executable). Note that both of these methods will ultimately call the command line executable. The advantage of the MATLAB script is that you can debug or plug in your own sorting methods. The advantage of the command line version is that it is independent of MATLAB and may be deployed on any UNIX system.

* Running from MATLAB:
Use ds001_sort_all.m

* Running from command line (preferred):
Be sure that the mountainsort/bin/mountainsort is in your PATH
> ./sort_tetrode1.sh
> ./sort_tetrode2.sh
> ./sort_ms11d45.sh
> ./sort_synth_EJ_K7.sh

