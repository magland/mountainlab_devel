function test_mec_01_29_2016_extract_raw

%extract_tetrode_data;
extract_probe_data

end

function extract_probe_data
mfile_path=fileparts(mfilename('fullpath'));
path_raw=sprintf('%s/../raw/entorhinal_ctx/probe',mfile_path);

raw_mat_fname=sprintf('%s/CH3_d04_e02_MECdata_new.mat',path_raw);
raw_mda_fname=sprintf('%s/CH3_d04_e02_MECdata_new.mda',path_raw);
probe1_fname=sprintf('%s/probe1_raw.mda',mfile_path);

if (~exist(raw_mda_fname,'file'))
    fprintf('Loading raw data...\n');
    L=load(raw_mat_fname);
    raw=L.MECdata.channelData';
    fprintf('Writing raw data...\n');
    writemda(raw,raw_mda_fname);
    writemda(L.MECdata.channelsMap,sprintf('%s/channels_map.mda',path_raw));
    writemda(L.MECdata.channelsBad,sprintf('%s/channels_bad.mda',path_raw));
end;

fprintf('Reading raw data...\n');
raw=readmda(raw_mda_fname);
channels_map=readmda(sprintf('%s/channels_map.mda',path_raw));
channels_bad=readmda(sprintf('%s/channels_bad.mda',path_raw));

channels_map'

%PROBE 1
fprintf('Writing probe1 data...\n');
probe=raw([1,2:19],1e6+1:35e6);
probe=probe(2:end,:)-repmat(probe(1,:),size(probe,1)-1,1);
writemda(probe,probe1_fname);

fprintf('Writing locations...\n');
L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];
writemda(L,[mfile_path,'/locations.mda']);


end

function extract_tetrode_data

mfile_path=fileparts(mfilename('fullpath'));
path_raw=sprintf('%s/../raw/entorhinal_ctx/tetrode',mfile_path);

raw_mat_fname=sprintf('%s/D09_d06_e04_MECdata_new.mat',path_raw);
raw_mda_fname=sprintf('%s/D09_d06_e04_MECdata_new.mda',path_raw);
tetrode3_fname=sprintf('%s/tetrode3_raw.mda',mfile_path);
tetrode4_fname=sprintf('%s/tetrode4_raw.mda',mfile_path);

if (~exist(raw_mda_fname,'file'))
    fprintf('Loading raw data...\n');
    L=load(raw_mat_fname);
    raw=L.MECdata.channelData';
    fprintf('Writing raw data...\n');
    writemda(raw,raw_mda_fname);
    writemda(L.MECdata.channelsMap,sprintf('%s/channels_map.mda',path_raw));
    writemda(L.MECdata.channelsBad,sprintf('%s/channels_bad.mda',path_raw));
end;

fprintf('Reading raw data...\n');
raw=readmda(raw_mda_fname);
channels_map=readmda(sprintf('%s/channels_map.mda',path_raw));
channels_bad=readmda(sprintf('%s/channels_bad.mda',path_raw));

channels_map'

%TETRODE 3
fprintf('Writing tetrode3 data...\n');
tetrode=raw([1,2:5],1e6+1:35e6);
tetrode=tetrode(2:end,:)-repmat(tetrode(1,:),size(tetrode,1)-1,1);
writemda(tetrode,tetrode3_fname);

%TETRODE 4
fprintf('Writing tetrode4 data...\n');
tetrode=raw([1,6:9],1e6+1:35e6);
tetrode=tetrode(2:end,:)-repmat(tetrode(1,:),size(tetrode,1)-1,1);
writemda(tetrode,tetrode4_fname);

fprintf('Writing locations...\n');
L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];
writemda(L,[mfile_path,'/locations.mda']);

end