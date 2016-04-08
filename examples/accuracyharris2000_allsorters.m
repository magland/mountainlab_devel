% accuracy vs IC on Harris 2000
% Barnett 4/8/16

clear; addpath('sorting_algs/demo_sort_001');
d = grab_harris2000_dataset;
o.detect_threshold = 3.5;   % sorter opts: threshold in stddev
o.detect_polarity = 'm';
o_acc.verb = 2; o_acc.usepre = 1;

sorters = {@simplesorter, @ds001_sort, @franklab_sort_2016_03_17_msdet4};
snames = {'simple','demosort001','franklab 3-17 det4'};
for i=[1 3]  %1:numel(sorters)
  fprintf('RUNNING sorter %s...\n',snames{i})
  [fk Q] = accuracy_anysorter_groundtrutheddata(sorters{i},d,o_acc,o);
  fprintf('Harris accuracy of sorter %s done; fk accuracies vs label k are... \n',snames{i})
  fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');
  fprintf('extended confusion matrix:\n')
  for i=1:size(Q,1), fprintf('%d\t',Q(i,:)), fprintf('\n\n'); end
end
