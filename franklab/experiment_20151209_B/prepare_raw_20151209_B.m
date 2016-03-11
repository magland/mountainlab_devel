function prepare_raw_20151209

mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode/20151209/20151209_s1_r1_s2_r2_s3_r3_tet9_ref.mat'];
L=load(raw_path);

list={'r1','r2','r3','s1','s2','s3'};

for ii=1:length(list)
    fprintf('ii=%d/%d\n',ii,length(list));
    str=list{ii};
    path0=[mfile_path,sprintf('/output_%s',str)];
    if (~exist(path0,'dir')) mkdir(path0); end;
    tmp=L.(str).channelData;
    tmp=tmp(:,1:4)-repmat(tmp(:,5),1,4);
    tmp=tmp';
    fprintf('Writing...\n');
    writemda(tmp,[path0,'/pre0.mda']);
end;
fprintf('Done.\n');

end