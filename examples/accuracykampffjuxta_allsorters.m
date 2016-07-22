% accuracy vs juxtacellular Kampff 32-ch MEA, 5-10 mins long
% Barnett 4/27/16, 7/22/16

clear; d = grab_kampffjuxta_dataset(1);  % choose Neto dataset 1 or 2

o_acc.verb = 3;  % 3 calls MV
o_acc.xc = 0; o_acc.usepre = 1;  % accuracy-testing opts

% common options for sorters
oc.detect_threshold = 2.5;   % sorter opts: threshold in stddev
oc.detect_polarity = 'm'; oc.sign=-1;
oc.elecadjmat = d.elecadjmat;    % branchv2 needs

sorters = {@simplesorter, @jfm_april_sort,@alg_scda_005};
snames = {'simple','jfm april','algscda005'};
for i=1:numel(sorters), os{i} = oc; end    % use common opts
% ... put any variations from common opts here:...
os{1}.num_fea=50;
os{2}.detectability_threshold=0;
os{3}.parfile = '/mnt/xfs1/home/ahb/ss_datasets/Kampff/kampff_M32.par'; % overrides
os{3}.eleccoordsfile = d.eleccoords;

testlist = [3]; %1:numel(sorters);   % which sorters to run

for i=testlist
  fprintf('RUNNING sorter %s...\n',snames{i})
  d.outdir = sprintf('%s/output%d',tempdir,i);  % separate output dirs for mv
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorters{i},d,o_acc,os{i});
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  if o_acc.verb, set(gcf,'name',snames{i}); end % label fig window
  fprintf('Kampff-accuracy of sorter %s done: fk accuracies vs label k are... \n',snames{i})
  fprintf('k   :'); fprintf('\t%d',1:numel(fk{i})), fprintf('\n');
  fprintf('f_k :'); fprintf('\t%.3f',fk{i}), fprintf('\n');
  fprintf('extended confusion matrix:\n')
  for j=1:size(Q{i},1), fprintf('%d\t',Q{i}(j,:)), fprintf('\n'); end
  fprintf('\n')
  fprintf('assuming one gndtruth neuron:   # correct %d, # missed %d, # false pos %d\n',Q{i}(1,1),sum(Q{i}(1,2:end)), Q{i}(2,1))
end

disp('ALL SORTERS KAMPFF-ACCURACY SUMMARY (SINGLE JC LABEL):')
for i=testlist, fprintf('sorter %s:\n',snames{i})
  fprintf('k   :'); fprintf('\t%d',1:numel(fk{i})), fprintf('\n');
  fprintf('f_k :'); fprintf('\t%.3f',fk{i}), fprintf('\n');
  fprintf(' assuming one gndtruth neuron:   # correct %d, # missed %d, # false pos %d\n',Q{i}(1,1),sum(Q{i}(1,2:end)), Q{i}(2,1))
end
