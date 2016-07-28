% accuracy vs IC on Martinez 2009 synth, for various sorters.
% Barnett 4/8/16, 4/19/16, algscda005 7/22/16 unfinished

clear; addpath('sorting_algs/demo_sort_001');
n = 4; d = grab_martinez2009_dataset(n);      % choose n=1..5 (see Martinez '09)
% (n=3 is most challenging, small ampl).
% Compare accuracies to Table 2 in their paper (given as false pos, missed).

o_acc.verb = 1; o_acc.xc = 1; o_acc.usepre = 1;  % accuracy-testing opts

% common options for sorters
oc.detect_threshold = 3.0;   % sorter opts: threshold in stddev

sorters = {@simplesorter, @ds001_sort, @franklab_sort_2016_03_17_msdet4, @jfm_april_sort, @validspike_wrapper, @alg_scda_005};
snames = {'simple','demosort001','franklab 3-17 det4','jfm april','validspike','algscda005'};
for i=1:numel(sorters), os{i} = oc; end    % use common opts
% ... put any variations from common opts here: ...
os{1}.detect_polarity = 'p';
% ds001 has no sign option
os{4}.sign = +1; os{4}.detectability_threshold = 0;  % no filtering of det score
os{4}.whiten=0;   % plain normalization (whitening breaks for M=1)
%os{5}.num_fea = 5;  % doesn't help much, still bad

testlist = [4]; %numel(sorters);   % which sorters to run (5 needs validspike)

for i=testlist
  fprintf('RUNNING sorter %s...\n',snames{i})
  d.outdir = sprintf('%s/output%d',tempdir,i);  % separate output dirs for mv
  [fk{i} Q{i}, perm{i}, info{i}] = accuracy_anysorter_groundtrutheddata(sorters{i},d,o_acc,os{i});
  if o_acc.verb, set(gcf,'name',snames{i}); end % label this fig window
  fprintf('Martinez-accuracy of sorter %s done: fk accuracies vs label k are... \n',snames{i})
  fprintf('k   :'); fprintf('\t%d',1:numel(fk{i})), fprintf('\n');
  fprintf('f_k :'); fprintf('\t%.3f',fk{i}), fprintf('\n');
  fprintf('extended confusion matrix:\n')
  for j=1:size(Q{i},1), fprintf('%d\t',Q{i}(j,:)), fprintf('\n'); end
  fprintf('\n')
end

disp('ALL SORTERS MARTINEZ-ACCURACY SUMMARY (SINGLE IC LABEL):')
for i=testlist, fprintf('sorter %s:\n',snames{i})
  fprintf('k   :'); fprintf('\t%d',1:numel(fk{i})), fprintf('\n');
  fprintf('f_k :'); fprintf('\t%.3f',fk{i}), fprintf('\n');
end

%diagnose_merge(info{4}.mergeinfo); % NB perm{i} and info{i} useful for diagnosing ith sorter
