Demo scripts for sorting and viewing output
Author: Jeremy Magland
Date: 4 Mar 2016

=== Step 1: Prepare raw data ===

You need to place the raw data in the appropriate locations within franklab/raw

Then run ds002_prepare_raw_data.m in MATLAB. This will create some output directories and place a pre0.mda file in each.


=== Step 2: Run the sorting ===

* Running from command line:
Be sure that the mountainsort/bin/mountainsort is in your PATH
> ./sort_tetrode1.sh
> ./sort_tetrode2.sh
> ./sort_ms11d45.sh
> ./sort_synth_EJ_K7.sh

