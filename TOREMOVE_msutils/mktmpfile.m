function fname=mktmpfile(X)

if (isstr(X)) fname=X; return; end; %already a file path! don't need to create a new file

mfile_path=fileparts(mfilename('fullpath'));
tmp_path=[mfile_path,'/../tmp'];

if ~exist(tmp_path,'dir') mkdir(tmp_path); end;

fname=sprintf('%s/tmp_%d_%d.mda',tmp_path,randi(100000),randi(100000));
writemda(X,fname);

end