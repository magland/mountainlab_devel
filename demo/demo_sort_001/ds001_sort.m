function [firings_path,pre_path]=ds001_sort(raw_path,output_path,sort_opts)
%DS001_SORT - Demo sorting script
%
% Syntax:  firings_path=ds001_sort(raw_path,output_path,sort_opts)
%
% Inputs:
%    raw_path - path to .mda of MxN raw signal data
%    output_path - path to existing DIRECTORY where all output will be written
%    sort_opts - (optional) sorting options, see def_sort_opts in this
%                   script
%
% Outputs:
%    firings_path - path to the firings.mda output file
%    pre_path - path to the preprocessed raw data file
%
% Other m-files required: isosplit2, mscmd_*

% Author: Jeremy Magland
% Mar 2016; Last revision: 1-Mar-2016

def_sort_opts.clip_size=100;
def_sort_opts.samplefreq=30000;
def_sort_opts.freq_min=300;
def_sort_opts.freq_max=6000;
def_sort_opts.outlier_threshold=500;
def_sort_opts.detect_threshold=3.5;
def_sort_opts.detect_interval=10;
def_sort_opts.individual_channels=1;
def_sort_opts.adjacency_matrix=[];
def_sort_opts.detectability_threshold=4;
def_sort_opts.shell_increment=0.5;
def_sort_opts.min_shell_size=100;

if nargin<3 sort_opts=struct; end;
sort_opts=ms_set_default_opts(sort_opts,def_sort_opts);

sort_opts.noise_subclusters.clip_size=sort_opts.clip_size;

for m=1:size(sort_opts.adjacency_matrix,1)
    sort_opts.adjacency_matrix(m,m)=0;
end;

path0=output_path;

o_filter.samplefreq=sort_opts.samplefreq;
o_filter.freq_min=sort_opts.freq_min;
o_filter.freq_max=sort_opts.freq_max;
o_filter.outlier_threshold=sort_opts.outlier_threshold;
o_detect.detect_threshold=sort_opts.detect_threshold;
o_detect.detect_interval=sort_opts.detect_interval;
o_detect.clip_size=sort_opts.clip_size;
o_detect.individual_channels=1;
o_detect.normalize=0;
o_noise_subclusters.detectability_threshold=sort_opts.detectability_threshold;
o_noise_subclusters.clip_size=sort_opts.clip_size;
o_branch_cluster.adjacency_matrix=sort_opts.adjacency_matrix;
o_branch_cluster.clip_size=sort_opts.clip_size;
o_branch_cluster.min_shell_size=sort_opts.min_shell_size;
o_branch_cluster.shell_increment=sort_opts.shell_increment;
o_branch_cluster.num_features=3;

%%%% Preprocessing
mscmd_bandpass_filter(raw_path,[path0,'/pre1.mda'],o_filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],struct);
mscmd_detect([path0,'/pre2.mda'],[path0,'/detect.mda'],o_detect);

%%%% Clustering
mscmd_branch_cluster_v1([path0,'/pre2.mda'],[path0,'/detect.mda'],'',[path0,'/firings1.mda'],o_branch_cluster);

%%%% Pruning
mscmd_remove_duplicates([path0,'/firings1.mda'],[path0,'/firings2.mda']);
mscmd_remove_noise_subclusters([path0,'/pre2.mda'],[path0,'/firings2.mda'],[path0,'/firings3.mda'],o_noise_subclusters);

%%%% Copying
mscmd_copy([path0,'/firings3.mda'],[path0,'/firings.mda']);

firings_path=[path0,'/firings.mda'];
pre_path=[path0,'/pre2.mda'];

end

