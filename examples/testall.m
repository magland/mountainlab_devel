% Script to run all unit tests and examples for ml_devel
% Barnett 4/6/16

clear; close all

% UNIT TESTS

loaddemowaveforms;
demo_dataset;
randomfirings
test_filter
test_detect_accuracy
detect_tests   % jfm

% SORTING ALGS

simplesorter

% VALIDATION

synthesize_timeseries
bestcolpermconfmat
times_labels_accuracy
times_labels_confusion_matrix
compare_two_sortings
accuracy_simplesorter
accuracy_ds001
accuracy_franksort_2016_03_17

% EXAMPLES WITH DEMO DATA

makesynthtimeseries
%driver_simplesorter   % already covered by simplesorter self-test

% EXAMPLES WITH NON-DEMO DATASETS

grab_EJnbhd_dataset
driver_sort_EJnbhd
grab_buzsaki_dataset
driver_sort_buzsaki
grab_harris2000_dataset

% ML PROCESSING

ms_detect3
mscmd_detect3

% THINGS IN DEVEL (SCRATCH) TO GO TO ML
addpath scratch_ahb

ms_detect4

if 0  % OBSOLETE THINGS STILL IN VARIOUS TESTERS
  ms_detect
end
