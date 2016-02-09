function experiment1_sort

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Raw extract options
opts.o_extract.num_channels=72;
opts.o_extract.channels=[37:52,68,69];
opts.o_extract.t1=0; opts.o_extract.t2=19e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial bandpass filter options
opts.o_filter.samplefreq=30000;
opts.o_filter.freq_min=600;
opts.o_filter.freq_max=4000;
opts.o_filter.outlier_threshold=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Whitening options
opts.o_whiten.ncomp=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-whitening bandpass filter options
opts.o_filter2.samplefreq=30000;
opts.o_filter2.freq_min=600;
opts.o_filter2.freq_max=4000;
opts.o_filter2.outlier_threshold=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detection options
opts.o_detect.inner_window_width=15;
opts.o_detect.outer_window_width=100000;
opts.o_detect.threshold=5;
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

mfile_path=fileparts(mfilename('fullpath'));
addpath(sprintf('%s/..',mfile_path));
outputdir_path=sprintf('%s/output',mfile_path);
if (~exist(outputdir_path,'dir')) mkdir(outputdir_path); end;
datfile_path=sprintf('%s/../raw/ms11d45.dat',mfile_path);
rawfile_path=sprintf('%s/output/raw.mda',mfile_path);

mscmd_extract(datfile_path,[outputdir_path,'/raw.mda'],opts.o_extract);

locations=get_frank_lab_locations;
opts.adjacency=ms_adjacency_matrix(locations,2);
writemda(locations,[outputdir_path,'/locations.mda']);

mountainsort_v1(rawfile_path,outputdir_path,opts);

end
