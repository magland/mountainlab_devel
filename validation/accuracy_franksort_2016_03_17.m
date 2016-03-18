% Showing how to use the accuracy-measuring function on a recent franklab expt
% sorter of jfm, but demo data.
% Barnett 3/18/16

clear
o = [];
o.detect_threshold = 3.0;   % ds001 uses thresh in sigma units
fk = accuracy_anysorter_groundtrutheddata(@franklab_sort_2016_03_17,[],[],o);
fprintf('franklab sort done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');

