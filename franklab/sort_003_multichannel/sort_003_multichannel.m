function [firings_path,pre_path]=sort_003_multichannel(raw_path,output_path,sort_opts)
%SORT_003_MULTICHANNEL - Version 003 of sorting based on shell method and isosplit2
%
% Syntax:  firings_path=sort_003_multichannel(raw_path,output_path,sort_opts)
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
% Other m-files required: isosplit2, mscmd_*, ms_*

% Author: Jeremy Magland
% Feb 2016; Last revision: 1-Mar-2016

if nargin<1 test_sort_003_multichannel; return; end;

def_sort_opts.clip_size=100;
def_sort_opts.filter.samplefreq=30000;
def_sort_opts.filter.freq_min=300;
def_sort_opts.filter.freq_max=6000;
def_sort_opts.filter.outlier_threshold=500;
def_sort_opts.detect.detect_threshold=3.5;
def_sort_opts.detect.detect_interval=10;
def_sort_opts.detect.individual_channels=1;
%def_sort_opts.min_cluster_size=10;
%def_sort_opts.detect.individual_channels=1;
def_sort_opts.adjacency_matrix=[];
def_sort_opts.noise_subclusters.detectability_threshold=4;
def_sort_opts.noise_subclusters.shell_increment=0.5;
def_sort_opts.noise_subclusters.min_shell_size=100;

if nargin<3 sort_opts=struct; end;
sort_opts=ms_set_default_opts(sort_opts,def_sort_opts);

sort_opts.noise_subclusters.clip_size=sort_opts.clip_size;

if (~sort_opts.detect.individual_channels)
    error('individual_channels must be set to 1');
end;

for m=1:size(sort_opts.adjacency_matrix,1)
    sort_opts.adjacency_matrix(m,m)=0;
end;

path0=output_path;

%%%% Preprocessing
mscmd_bandpass_filter(raw_path,[path0,'/pre1.mda'],sort_opts.filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],struct);
mscmd_detect([path0,'/pre2.mda'],[path0,'/detect.mda'],sort_opts.detect);

%%%% Clustering
mscmd_branch_cluster_v1([path0,'/pre2.mda'],[path0,'/detect.mda'],'',[path0,'/firings1.mda']);

%%%% Pruning
mscmd_remove_duplicates([path0,'/firings1.mda'],[path0,'/firings2.mda']);
mscmd_remove_noise_subclusters([path0,'/pre2.mda'],[path0,'/firings2.mda'],[path0,'/firings3.mda'],sort_opts.noise_subclusters);

%%%% Copying
mscmd_copy([path0,'/firings3.mda'],[path0,'/firings.mda']);

firings_path=[path0,'/firings.mda'];
pre_path=[path0,'/pre2.mda'];

end

function test_extract_raw_data(raw_path,output_path,tetrode_num)

raw_mat_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mat',raw_path);
raw_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mda',raw_path);
tetrode_fname=sprintf('%s/pre0.mda',output_path);

if (~exist(raw_mda_fname,'file'))
    fprintf('Loading raw data...\n');
    L=load(raw_mat_fname);
    raw=L.dl12_20151208_NNF_r1_tet16_17.channelData';
    fprintf('Writing raw data...\n');
    writemda(raw,raw_mda_fname);
end;

if (~exist(tetrode_fname,'file'))
    fprintf('Reading raw data...\n');
    raw=readmda(raw_mda_fname);

    if (tetrode_num==1)
        tetrode=raw([1,3:6],(1e6+1):26e6);
        tetrode=tetrode(2:end,:)-repmat(tetrode(1,:),size(tetrode,1)-1,1);
    elseif (tetrode_num==2)
        tetrode=raw([1,7:10],(1e6+1):26e6);
        tetrode=tetrode(2:end,:)-repmat(tetrode(1,:),size(tetrode,1)-1,1);
    end;
    fprintf('Writing tetrode data...\n');
    writemda(tetrode,tetrode_fname);
end;

L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];
writemda(L,[output_path,'/locations.mda']);

end

function test_sort_003_multichannel

close all;

tetrode_num=1;

%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode'];
path0=[mfile_path,sprintf('/output_tetrode%d',tetrode_num)];
if ~exist(path0,'dir') mkdir(path0); end;

%%%% Extract raw data
test_extract_raw_data(raw_path,path0,tetrode_num);

%%%% Sort
[firings_path,pre_path]=sort_003_multichannel([path0,'/pre0.mda'],path0);

%%%% View output
mv.mode='overview2';
mv.raw=pre_path;
mv.firings=firings_path'
ms_mountainview(mv);

end


