function test_hippocampal_01_28_2016_extract_raw

mfile_path=fileparts(mfilename('fullpath'));
path_raw=sprintf('%s/../raw/hippocampal/tetrode',mfile_path);
path0=sprintf('%s/output',mfile_path);
if ~exist(path0,'dir')
    mkdir(path0);
end;

raw_mat_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mat',path_raw);
raw_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mda',path_raw);
tetrode1_fname=sprintf('%s/tetrode1_raw.mda',mfile_path);
tetrode2_fname=sprintf('%s/tetrode2_raw.mda',mfile_path);

if (~exist(raw_mda_fname,'file'))
    fprintf('Loading raw data...\n');
    L=load(raw_mat_fname);
    raw=L.dl12_20151208_NNF_r1_tet16_17.channelData';
    fprintf('Writing raw data...\n');
    writemda(raw,raw_mda_fname);
end;

fprintf('Reading raw data...\n');
raw=readmda(raw_mda_fname);

tetrode1=raw([1,3:6],(1e6+1):26e6);
tetrode1=tetrode1(2:end,:)-repmat(tetrode1(1,:),size(tetrode1,1)-1,1);
fprintf('Writing tetrode1 data...\n');
writemda(tetrode1,tetrode1_fname);

tetrode2=raw([1,7:10],(1e6+1):26e6);
tetrode2=tetrode2(2:end,:)-repmat(tetrode2(1,:),size(tetrode2,1)-1,1);
fprintf('Writing tetrode2 data...\n');
writemda(tetrode2,tetrode2_fname);

L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];
writemda(L,[mfile_path,'/locations.mda']);


end