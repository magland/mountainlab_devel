% Showing how to use the accuracy-measuring function on a recent franklab expt
% sorter of jfm, but demo data. Note all opts are built into the sorter.
% Barnett 3/18/16, 4/7/16

clear
%fk = accuracy_anysorter_groundtrutheddata(@franklab_sort_2016_03_17,[],[]);
fk = accuracy_anysorter_groundtrutheddata(@franklab_sort_2016_03_17_msdet4,[],[]);
fprintf('franklab sort done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');

