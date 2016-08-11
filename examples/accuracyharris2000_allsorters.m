% accuracy vs IC on Harris 2000, for various sorters
% Barnett 4/8/16

clear; addpath('sorting_algs/demo_sort_001');
d = grab_harris2000_dataset;

o_acc.verb = 0; o_acc.xc = 0; o_acc.usepre = 1;  % accuracy-testing opts

% common options for sorters
oc.detect_threshold = 3.0;   % sorter opts: threshold in stddev
oc.detect_polarity = 'm'; oc.sign=-1;

sorters = {@simplesorter, @ds001_sort, @franklab_sort_2016_03_17_msdet4, @jfm_april_sort, @validspike_wrapper, @alg_scda_005};
snames = {'simple','demosort001','franklab 3-17 det4','jfm april','validspike','algscda005'};
for i=1:numel(sorters), os{i} = oc; end    % use common opts
% ... put any variations from common opts here:...
%os{2}.detect_threshold = 3.0;
%os{3}.detect_threshold = 3.0;
os{4}.detectability_threshold=5;

os{6}.detect_threshold=3.5; os{6}.clip_size=20; os{6}.use_whitening=0; os{6}.use_mask_out_artifacts=0; os{6}.freq_max=inf; os{6}.num_fea=20; %os{6}.neglogprior=30;


testlist = 6; %[1 4 5 6]; %1:4; %numel(sorters);   % which sorters to run

for i=testlist
  fprintf('RUNNING sorter %s...\n',snames{i})
  d.outdir = sprintf('%s/output%d',tempdir,i);  % separate output dirs for mv
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorters{i},d,o_acc,os{i});
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
  % miss and fp rates for the case of 1 gnd-truth...
  mk{i} = sum(Q{i}(1,2:end))/Q{i}(1,1);   % frac missed (<=1)
  pk{i} = Q{i}(2,1)/Q{i}(1,1);          % false pos ratio to good
  ek{i} = Q{i}(2,1)/sum(Q{i}(:,1));   % est false pos frac (as in Ekanadham, <=1)
  fprintf('misfrac'); fprintf('\t%.3f',mk{i}), fprintf('\n');
  fprintf('fp/good'); fprintf('\t%.3f',pk{i}), fprintf('\n');
  fprintf('estfpfr'); fprintf('\t%.3f',ek{i}), fprintf('\n');
end

%diagnose_merge(info{4}.mergeinfo); % NB perm{i} and info{i} useful for diagnosing ith sorter. mergeinfo only exists for old sorters.

% show the ith sorting...
i=6; outdir = sprintf('%s/output%d',tempdir,i);
if 1
mv.mode='overview2'; mv.raw = d.timeseries;   % in MV...
mv.filt=info{i}.filtfile; mv.pre=info{i}.prefile;mv.samplerate=d.samplerate;
mv.firings= [outdir,'/firings_permed.mda'];mountainview(mv);

Y=readmda([outdir,'/pre2.mda']); F=readmda([outdir,'/firings_permed.mda']);
Ft=readmda(d.truefirings); spikespy({Y,F(2,:),F(3,:),'sorted'},{Y,Ft(2,:),Ft(3,:),'truth'});
end

if 0  % show IC signal too just to check gnd truth worked.
  YIC=readmda('/tmp/harrisoutput/harris2000_IC.mda');
Y=readmda([outdir,'/pre2.mda']); F=readmda([outdir,'/firings_permed.mda']);
Ft=readmda(d.truefirings); spikespy({Y,F(2,:),F(3,:),'sorted'},{Y,Ft(2,:),Ft(3,:),'truth'},{YIC,Ft(2,:),Ft(3,:),'IC'});
end