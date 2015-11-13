function ms_example1_process

addpath([fileparts(mfilename('fullpath')),'/..']);

opts.locations=get_frank_lab_locations;

%input_file_path='dev/test/testA.dat';
%output_dir_path='test_output';
input_file_path=[fileparts(mfilename('fullpath')),'/../example1_data/ms11d45.dat'];
%output_dir_path=[fileparts(mfilename('fullpath')),'/../example1_output'];
output_dir_path=[fileparts(mfilename('fullpath')),'/../example1c_output'];

opts.num_channels=72;
opts.channels=[37:52,68,69];
%opts.timepoints=1:1e6;
opts.timepoints=[1:2.9e6,3.0e6:19e6]; %Something seems to change around timepoint 19-20 million
opts.dtype='int16';
opts.samplefreq=30000;
opts.freq_min=300;
opts.freq_max=10000;

%comment this out if you change any of the above parameters that
%relate to the pre-processing so the raw.mda file can be re-generated
opts.preprocessed_input_file_path=[fileparts(mfilename('fullpath')),'/../example1_data/ms11d45A_pre.mda'];

opts.detect_interval=40;
opts.detect_threshold=5;
opts.clip_size=80;
opts.num_pca_features=6;
opts.adjacency_radius=2;
opts.min_cluster_size=20;
opts.cross_correlogram_max_dt=1500;
opts.whiten=1;
opts.num_whitening_components=3;
opts.cluster_outlier_alpha=0.01;

mountainsort(input_file_path,output_dir_path,opts);

%ms_example1_view

addpath([fileparts(mfilename('fullpath')),'/..']);
addpath([fileparts(mfilename('fullpath')),'/../view']);
mountainview(output_dir_path);
view_cross_correlograms([output_dir_path,'/cross-correlograms.mda'],0);

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
