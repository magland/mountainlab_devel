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
    raw=raw(:,(1e6+1):26e6);
    fprintf('Writing raw subset data...\n');
    writemda(raw,raw_subset_mda_fname);
end;

opts_pre=get_sorting_options;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial bandpass filter options
opts_pre.o_filter.samplefreq=30000;
opts_pre.o_filter.freq_min=600;
opts_pre.o_filter.freq_max=4000;
opts_pre.o_filter.outlier_threshold=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Whitening options
opts_pre.o_whiten.ncomp=2;
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
mountainsort_v1('',path0,opts_sort);

end

function opts=get_sorting_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detection options
opts.o_detect.inner_window_width=15;
opts.o_detect.outer_window_width=1000;
opts.o_detect.threshold=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature extraction options
opts.o_features.num_features=6;
opts.o_features.clip_size=60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clustering options
opts.o_cluster=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster splitting options
opts.o_split_clusters.num_features=3; 
opts.o_split_clusters.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Template extraction options
opts.o_templates.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster consolidation options
opts.o_consolidate.coincidence_threshold=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit options
opts.o_fit=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross-correlogram options
opts.o_cross_correlograms.max_dt=1500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clip extraction options
opts.o_extract_clips.clip_size=opts.o_templates.clip_size;
end