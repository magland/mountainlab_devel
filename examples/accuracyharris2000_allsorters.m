% accuracy vs IC on Harris 2000, for various sorters
% Barnett 4/8/16

clear; addpath('sorting_algs/demo_sort_001');
d = grab_harris2000_dataset;

o_acc.verb = 1; o_acc.usepre = 1;  % accuracy-testing opts

% common options for sorters
oc.detect_threshold = 3.0;   % sorter opts: threshold in stddev
oc.detect_polarity = 'm'; oc.sign=-1;

sorters = {@simplesorter, @ds001_sort, @franklab_sort_2016_03_17_msdet4, @jfm_april_sort, @validspike_wrapper};
snames = {'simple','demosort001','franklab 3-17 det4','jfm april','validspike'};
for i=1:numel(sorters), os{i} = oc; end    % use common opts
% ... put any variations from common opts here:...
%os{2}.detect_threshold = 3.0;
%os{3}.detect_threshold = 3.0;

testlist = 1:4; %numel(sorters);   % which sorters to run

for i=testlist
  fprintf('RUNNING sorter %s...\n',snames{i})
  d.outdir = sprintf('%s/output%d',tempdir,i);  % separate output dirs for mv
  [fk{i} Q{i}] = accuracy_anysorter_groundtrutheddata(sorters{i},d,o_acc,os{i});
  if o_acc.verb, set(gcf,'name',snames{i}); end % label fig window
  fprintf('Harris-accuracy of sorter %s done: fk accuracies vs label k are... \n',snames{i})
  fprintf('k   :'); fprintf('\t%d',1:numel(fk{i})), fprintf('\n');
  fprintf('f_k :'); fprintf('\t%.3f',fk{i}), fprintf('\n');
  fprintf('extended confusion matrix:\n')
  for j=1:size(Q{i},1), fprintf('%d\t',Q{i}(j,:)), fprintf('\n'); end
  fprintf('\n')
end

disp('ALL SORTERS HARRIS-ACCURACY SUMMARY (SINGLE IC LABEL):')
for i=testlist, fprintf('sorter %s:\n',snames{i})
  fprintf('k   :'); fprintf('\t%d',1:numel(fk{i})), fprintf('\n');
  fprintf('f_k :'); fprintf('\t%.3f',fk{i}), fprintf('\n');
end
