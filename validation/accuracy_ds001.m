% Showing how to use the accuracy-measuring function on demo sorter 001
% Barnett 3/17/16

clear; addpath('sorting_algs/demo_sort_001');
o = [];
o.detect_threshold = 3.0;   % ds001 uses thresh in sigma units
fk = accuracy_anysorter_groundtrutheddata(@ds001_sort,[],[],o);
fprintf('demo_sort_001 done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');
