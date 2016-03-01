function ds001_prepare_raw_data

prepare_hippocampal_tetrode_data(1);
prepare_hippocampal_tetrode_data(2);

end

function prepare_hippocampal_tetrode_data(tetrode_num)

fprintf('prepare_hippocampal_tetrode_data(%d)\n',tetrode_num);

mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../../franklab/raw/hippocampal/tetrode'];
output_path=[mfile_path,sprintf('/output_tetrode%d',tetrode_num)];
if ~exist(output_path,'dir') mkdir(output_path); end;

raw_mat_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mat',raw_path);
raw_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mda',raw_path);
tetrode_fname=sprintf('%s/pre0.mda',output_path);

fprintf('Loading raw data...\n');
L=load(raw_mat_fname);
raw=L.dl12_20151208_NNF_r1_tet16_17.channelData';
fprintf('Writing raw data...\n');
writemda(raw,raw_mda_fname);

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

L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];
writemda(L,[output_path,'/locations.mda']);

end
