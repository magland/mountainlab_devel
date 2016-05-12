% accuracy vs juxtacellular Kampff 32-ch MEA, 5-10 mins long
% Barnett 4/27/16

clear; d = grab_kampffjuxta_dataset(2);

o_acc.verb = 3; o_acc.xc = 0; o_acc.usepre = 1;  % accuracy-testing opts

% common options for sorters
oc.detect_threshold = 2.5;   % sorter opts: threshold in stddev
oc.detect_polarity = 'm'; oc.sign=-1;
oc.elecadjmat = d.elecadjmat;    % branchv2 needs

sorters = {@simplesorter, @jfm_april_sort};
snames = {'simple','jfm april'};
for i=1:numel(sorters), os{i} = oc; end    % use common opts
% ... put any variations from common opts here:...
os{1}.num_fea=50;
os{2}.detectability_threshold=0;

testlist = 2; %1:numel(sorters);   % which sorters to run

for i=testlist
  fprintf('RUNNING sorter %s...\n',snames{i})
  d.outdir = sprintf('%s/output%d',tempdir,i);  % separate output dirs for mv
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorters{i},d,o_acc,os{i});
  if o_acc.verb, set(gcf,'name',snames{i}); end % label fig window
  fprintf('Kampff-accuracy of sorter %s done: fk accuracies vs label k are... \n',snames{i})
  fprintf('k   :'); fprintf('\t%d',1:numel(fk{i})), fprintf('\n');
  fprintf('f_k :'); fprintf('\t%.3f',fk{i}), fprintf('\n');
  fprintf('extended confusion matrix:\n')
  for j=1:size(Q{i},1), fprintf('%d\t',Q{i}(j,:)), fprintf('\n'); end
  fprintf('\n')
end

disp('ALL SORTERS KAMPFF-ACCURACY SUMMARY (SINGLE JC LABEL):')
for i=testlist, fprintf('sorter %s:\n',snames{i})
  fprintf('k   :'); fprintf('\t%d',1:numel(fk{i})), fprintf('\n');
  fprintf('f_k :'); fprintf('\t%.3f',fk{i}), fprintf('\n');
end
