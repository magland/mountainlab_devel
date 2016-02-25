% Basic demo of spike sorting the shipped demo data with mountainsort pipeline.
% Very alpha. *** indicates to fix.
% Barnett 2/19/16

if ~exist('unit_tests/demo_data/demotimeseries.mda'), writedemotimeseries; end
clear
% (no filtering since the demo waveforms were already filtered & zero-mean)

% detect events as times on each channel...
det_opts.threshold = 80;          % absolute magnitude
det_opts.individual_channels = 0;  % take abs max across all channels
path0 = 'unit_tests/demo_data';
raw = [path0 '/demotimeseries.mda'];
mscmd_detect(raw,[path0,'/detect.mda'],det_opts);
CT = readmda([path0,'/detect.mda']); times=CT(2,:); clear CT
fprintf('detect found %d events\n',numel(times))
%Y = readmda(raw); spikespy({Y,times,0*times,'Y & detect'});  % 0 = no labels

clip_opts.clip_size=40;    % 2 ms, enough
mscmd_extract_clips(raw,[path0,'/detect.mda'],[path0,'/clips.mda'],clip_opts);
%X = readmda([path0,'/clips.mda']); spikespy({X,'clips X'});

fprintf('read clips...\n'); clips=readmda([path0,'/clips.mda']);
% Cluster the clips... AHB attempt to make a simple script
[M,T,L]=size(clips);
clus_opts.num_features=12;
iso_opts.whiten_at_each_comparison = 0;   % iso2 fails if 1  ***
X = ms_event_features(clips,clus_opts.num_features);  % for viz or clustering
clus = 1;   % clustering style - todo: should be opts via a cluster func!
if clus==0  % plain isosplit on *raw* clips
  fullX = reshape(clips,[M*T,L]); [labels,info] = isosplit2(fullX,iso_opts); 
  peaks = 0*labels;    % dummy peak ampls
elseif clus==1   % dimension reduction via PCA, then isosplit
  [labels,info] = isosplit2(X,iso_opts); 
  peaks = 0*labels;    % dummy peak ampls
elseif clus==2   % .. or shell stuff.   Can still crash
  clus_opts.min_shell_count=2000;
  clus_opts.shell_increment=1;
  clus_opts.merge_threshold=0.8;
  fprintf('--- jfm shell cluster ---\n');
  [labels,peaks]=jfm_shell_cluster(clips,clus_opts);  % isosplit2 fails ***
end
% write out firings... (for MV and for conf-mat)
K=max(labels);
pops = histc(labels,1:K)
firings=zeros(4,L);
detect=readmda([path0,'/detect.mda']);
times=detect(2,:);
firings(1:2,:)=detect; firings(3,:)=labels; firings(4,:)=peaks;
writemda(firings,[path0,'/firings.mda']);

% view...
plot_labeled_pts(X,labels);  % visualize with feature space even if not used
Y = readmda(raw); spikespy({Y,times,labels,'Y sorted'});  % Ctrl-E fails ***
% how do I now best look at labeled clips? try MV  *** needs pops below W's
view_params.mode='overview2';
view_params.raw=raw;
view_params.firings=[path0,'/firings.mda'];
view_params.sampling_freq = 2e4;   % determined by default waveforms
ms_mountainview(view_params);    % *** no output if peaks=nan

fprintf('compute accuracy vs known ground-truth...\n')
max_matching_offset = 10;    % in samples; 0.5 ms
mscmd_confusion_matrix([path0 '/truefirings.mda'],[path0 '/firings.mda'],[path0 '/accconfmat.mda'],max_matching_offset);
Q = readmda([path0 '/accconfmat.mda'])     % is this the best-perm'ed Q? no ***
figure; imagesc(Q); colorbar; ylabel('truth'); xlabel('sorted');
title('extended accuracy confusion matrix');
% clus=0,1 give false splits for high-pop clusters; merge of lowest 2

% now generate accuracy summary...
