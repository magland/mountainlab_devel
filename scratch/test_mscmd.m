function test_mscmd

addpath('mountainview/src/spikespy/matlab');
addpath('processing');
addpath('util');

path0='example1_data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_extract.num_channels=72;
o_extract.channels=[37:52,68,69];
%o_extract.t1=0; o_extract.t2=1e6;
o_extract.t1=0; o_extract.t2=19e6;
%o_extract.t1=0; o_extract.t2=5e7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_filter.samplefreq=30000;
o_filter.freq_min=300;
o_filter.freq_max=4000;
o_filter.outlier_threshold=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_whiten.ncomp=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_filter2.samplefreq=30000;
o_filter2.freq_min=600;
o_filter2.freq_max=10000;
o_filter2.outlier_threshold=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_detect.inner_window_width=40;
o_detect.outer_window_width=100000;
o_detect.threshold=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_features.num_features=6;
o_features.clip_size=200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_cluster=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_split_clusters.num_features=3;
o_split_clusters.clip_size=o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_remove_outliers=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_templates.clip_size=o_features.clip_size;

% locations=get_frank_lab_locations;
% AM=ms_adjacency_matrix(locations,2);
% writemda(AM,[path0,'/adjacency.mda']);
% writemda(locations,[path0,'/locations.mda']);

timerA=tic;
mscmd_extract([path0,'/ms11d45.dat'],[path0,'/raw.mda'],o_extract);
mscmd_bandpass_filter([path0,'/raw.mda'],[path0,'/filt.mda'],o_filter);
mscmd_whiten([path0,'/filt.mda'],[path0,'/filt_white.mda'],o_whiten);
mscmd_bandpass_filter([path0,'/filt_white.mda'],[path0,'/filt2_white.mda'],o_filter2);
mscmd_detect([path0,'/filt2_white.mda'],[path0,'/detect.mda'],o_detect);
mscmd_features([path0,'/filt2_white.mda'],[path0,'/detect.mda'],[path0,'/adjacency.mda'],[path0,'/features.mda'],o_features);
mscmd_cluster([path0,'/features.mda'],[path0,'/cluster.mda'],o_cluster);
mscmd_templates([path0,'/filt2_white.mda'],[path0,'/cluster.mda'],[path0,'/templates.mda'],o_templates);
mscmd_consolidate([path0,'/cluster.mda'],[path0,'/templates.mda'],[path0,'/cluster0.mda'],[path0,'/templates0.mda'],[path0,'/load_channels0.mda']);
mscmd_templates([path0,'/raw.mda'],[path0,'/cluster0.mda'],[path0,'/templates0_raw.mda'],o_templates);

%mscmd_split_clusters([path0,'/filt2_white.mda'],[path0,'/cluster.mda'],[path0,'/cluster2.mda'],o_split_clusters);
%mscmd_remove_outliers([path0,'/filt2_white.mda'],[path0,'/cluster2.mda'],[path0,'/cluster3.mda'],o_remove_outliers);
%mscmd_templates([path0,'/filt2_white.mda'],[path0,'/cluster3.mda'],[path0,'/templates.mda'],o_templates);
%mscmd_consolidate([path0,'/cluster3.mda'],[path0,'/templates.mda'],[path0,'/cluster0.mda'],[path0,'/templates0.mda'],[path0,'/load_channels0.mda']);
%mscmd_templates([path0,'/raw.mda'],[path0,'/cluster0.mda'],[path0,'/templates_raw.mda'],o_templates);

fprintf('Elapsed time: %g sec',toc(timerA));

% ch=3;
% out_path=[path0,sprintf('/sort-%d.mda',ch)];
% mscmd_sort([path0,'/filt2_white.mda'],[path0,'/detect.mda'],out_path,ch,o_sort);

% X=readmda([path0,'/filt2_white.mda']);
% detect=readmda([path0,'/detect.mda']);
% T=detect(2,:);
% L=detect(1,:);
% spikespy({X,T,L});

%X=readmda_data_beginning([path0,'/test1.mda'],3e6);
%spikespy(X);

%X=readmda_data_beginning([path0,'/test3.mda'],3e6);
%spikespy(X);

end

function mscmd_sort(X_path,detect_path,out_path,ch,opts)

locations=get_frank_lab_locations;
locations
AM=ms_adjacency_matrix(locations,2);
figure; imagesc(AM);
AM(ch,ch)=0;
neighbor_channels=[ch,find(AM(ch,:))]
return;

disp('reading...');
X=readmda(X_path);
X=X(neighbor_channels,:);
detect=readmda(detect_path);

inds=find(detect(1,:)==ch);
times=detect(2,inds);
disp('extract clips...');
clips=ms_extract_clips(X,times,opts.clip_size);
disp('features...');
FF=ms_event_features(clips,opts.num_features);
disp('isosplit...');
labels=isosplit(FF);
disp('view...');
ss_view_clusters(FF,labels);

spikespy({X,times,labels});

end

function L=get_frank_lab_locations

L=[...
0.0,0.0;...
-0.5,1.0;...
0.5,1.0;...
-1.0,2.0;...
1.0,2.0;...
-1.0,3.0;...
1.0,3.0;...
-1.0,4.0;...
1.0,4.0;...
-1.0,5.0;...
1.0,5.0;...
-1.0,6.0;...
1.0,6.0;...
-1.0,7.0;...
1.0,7.0;...
-1.0,8.0;...
1.0,8.0;...
-1.0,9.0;...
];

end
