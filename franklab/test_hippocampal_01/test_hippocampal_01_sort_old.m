function test_hippocampal_01_sort

mfile_path=fileparts(mfilename('fullpath'));
path_raw=sprintf('%s/../raw/hippocampal/tetrode',mfile_path);
path0=sprintf('%s/output',mfile_path);
if ~exist(path0,'dir')
    mkdir(path0);
end;

raw_mat_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mat',path_raw);
raw_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mda',path_raw);
raw_subset_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17_subset.mda',path_raw);

if (~exist(raw_mda_fname,'file'))
    fprintf('Loading raw data...\n');
    L=load(raw_mat_fname);
    raw=L.dl12_20151208_NNF_r1_tet16_17.channelData';
    fprintf('Writing raw data...\n');
    writemda(raw,raw_mda_fname);
    %spikespy({raw});
    raw=raw(3:6,(1e6+1):26e6);
    fprintf('Writing raw subset data...\n');
    writemda(raw,raw_subset_mda_fname);
end;

writemda(get_geometry,sprintf('%s/locations.mda',path0));

opts_pre=get_sorting_options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial bandpass filter options
opts_pre.o_filter.samplefreq=30000;
opts_pre.o_filter.freq_min=600;
opts_pre.o_filter.freq_max=4000;
opts_pre.o_filter.outlier_threshold=500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Whitening options
opts_pre.o_whiten.ncomp=1; %Use 0 for second tetrode, 1 for first tetrode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-whitening bandpass filter options
opts_pre.o_filter2.samplefreq=30000;
opts_pre.o_filter2.freq_min=600;
opts_pre.o_filter2.freq_max=4000;
opts_pre.o_filter2.outlier_threshold=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mscmd_bandpass_filter(raw_subset_mda_fname,[path0,'/pre1.mda'],opts_pre.o_filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],opts_pre.o_whiten);
mscmd_bandpass_filter([path0,'/pre2.mda'],[path0,'/pre3.mda'],opts_pre.o_filter2);

opts_sort=get_sorting_options;
opts_sort.adjacency=ones(10,10);
tmp_mountainsort('',path0,opts_sort);

end

function opts=get_sorting_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detection options
opts.o_detect.inner_window_width=15;
opts.o_detect.outer_window_width=1000;
opts.o_detect.threshold=3;
opts.o_detect.use_pre1=0;
opts.o_detect.individual_channels=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature extraction options
opts.o_features.num_features=6;
opts.o_features.clip_size=60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clustering options
opts.o_cluster.ks_threshold=1.0;
opts.o_cluster.K_init=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster splitting options
% opts.o_split_clusters.num_features=3; 
% opts.o_split_clusters.clip_size=opts.o_features.clip_size;
% opts.o_split_clusters.use_pre1=1;
% opts.o_split_clusters.ks_threshold=opts.o_cluster.ks_threshold;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Template extraction options
opts.o_templates.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster consolidation options
%opts.o_consolidate.coincidence_threshold=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit options
opts.o_fit=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross-correlogram options
opts.o_cross_correlograms.max_dt=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clip extraction options
opts.o_extract_clips.clip_size=opts.o_templates.clip_size;
end

function L=get_geometry
L=[...
-1.0,5;...
1,5;...
-2,3;...
-2,2;...
-2,1;...
-2,0;...
2,3;...
2,2;...
2,1;...
2,0;...
];

end
