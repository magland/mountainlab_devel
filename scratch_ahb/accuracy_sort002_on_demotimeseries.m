function accuracy_sort002_on_demotimeseries
% measure accuracy vs known firing times on synthetic demo data.
% Barnett 2/25/16


%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'))

addpath([mfile_path,'/../franklab/sort_002_multichannel']);
dpath = [mfile_path,'/../unit_tests/demo_data'];
raw=[dpath,'/demotimeseries.mda'];
if ~exist(raw), writedemotimeseries; end

mfile_path=fileparts(mfilename('fullpath'));
output_path=[mfile_path,'/output'];
if ~exist(output_path,'dir') mkdir(output_path); end;

%%%% Sort
sort_opts=struct;
sort_opts.detect.detect_threshold=4;   % in sigmas; don't vary much
sort_opts.detect.detect_interval=15;
sort_opts.detectability_threshold=0; % controls exclusion of noise clusters
sort_opts.plausibility_threshold=inf;  % exc
sort_opts.clip_size=50;          % seems to have big effect? (truth at most 46)
[firings_path,pre_path]=sort_002_multichannel(raw,output_path,sort_opts);

if 0 %%%% View output
  mv.mode='overview2';
  mv.raw=pre_path;
  mv.firings=firings_path;
  ms_mountainview(mv);
end

firings=readmda(firings_path);
times=firings(2,:); labels=firings(3,:);
K = max(labels);  % # claimed neurons
ii = labels>0;
if 0                % kill unclassified
  times = times(ii); labels=labels(ii);
else                % treat unclass as a new type
  labels(~ii) = K+1; K = max(labels);
end
Y=readmda(pre_path);                       % how preprocessed is this?
clips=ms_extract_clips(Y,times,60);
templates=ms_templates(clips,labels);
%figure; ms_view_templates(templates); % jfm W view

% load the truth...
Wtrue = readmda([dpath,'/waveforms_EJ_2005-04-26_elec359_20kHz_K7.mda']);
plot_spike_shapes(Wtrue,'true W');
truefirings=readmda([dpath,'/truefirings.mda']);
truetimes=truefirings(2,:); truelabels=truefirings(3,:);

fprintf('compute accuracy vs known ground-truth...\n')
o.max_matching_offset = 10;    % in samples; 0.5 ms
if 1     % old matlab validspike way
  [perm Q acc t] = times_labels_accuracy(truetimes,truelabels,times,labels,o);
else     % C fast way, but not full times-labels best search
  mscmd_confusion_matrix([dpath,'/truefirings.mda'],firings_path,[output_path '/accconfmat.mda'],o.max_matching_offset);
  Q = readmda([output_path '/accconfmat.mda']);   % un-permed
  [perm Q] = bestcolpermconfmat(Q);
  disp('best permed extended confusion matrix:'), Q
  [~,acc] = labels_similarity(Q);
end
fk = acc.p;   % our accuracy measure

figure; imagesc(Q); colorbar; ylabel('truth'); xlabel('sorted');
title('best extended accuracy confusion matrix');

[~,iperm] = sort(perm(1:K));
plot_spike_shapes(templates(:,:,iperm),'best-permed W');

addpath ~/spikespy/matlab/  % prefer old spikespy
spikespy({Y,truetimes,truelabels,'Y, truth'},{Y,times,perm(labels),'Y, sorted'});

keyboard    % use dbcont to continue
end
