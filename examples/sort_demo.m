% Demo of spike sorting the demo data with mountainsort.
% Barnett 2/19/16

if ~exist('unit_tests/demo_data/demotimeseries.mda'), writedemotimeseries; end
clear
% no filtering since the demo waveforms were already filtered, zero-mean

% detect events as times on each channel...
det_opts.threshold = 100;          % absolute magnitude
det_opts.individual_channels = 0;  % take abs max across all channels
path0 = 'unit_tests/demo_data';
raw = [path0 '/demotimeseries.mda'];
mscmd_detect(raw,[path0,'/detect.mda'],det_opts);
CT = readmda([path0,'/detect.mda']); times=CT(2,:); clear CT
fprintf('detect found %d events\n',numel(times))
%Y = readmda(raw); spikespy({Y,times,0*times,'Y & detect'});  % 0 = no labels

clip_opts.clip_size=40;
mscmd_extract_clips(raw,[path0,'/detect.mda'],[path0,'/clips.mda'],clip_opts);
%X = readmda([path0,'/clips.mda']); spikespy({X,'clips X'});

fprintf('read clips...\n'); clips=readmda([path0,'/clips.mda']);
[M,T,L]=size(clips);

shell_opts.min_shell_count=2000;
shell_opts.shell_increment=1;
shell_opts.num_features=12;
shell_opts.merge_threshold=0.8;
fprintf('--- jfm shell cluster ---\n');
[labels,peaks]=jfm_shell_cluster(clips,shell_opts);  % isosplit2 fails ***
K=max(labels);
firings=zeros(4,L);
detect=readmda([path0,'/detect.mda']);
times=detect(2,:);
firings(1:2,:)=detect;
firings(3,:)=labels;
firings(4,:)=peaks;
writemda(firings,[path0,'/firings.mda']);

fprintf('computing accuracy vs known ground-truth...\n')
% mscmd_confusion_matrix w/ truefirings.mda vs firings.mda
% ...
