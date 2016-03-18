function prepare_epochs_s1r1

mfile_path=fileparts(mfilename('fullpath'));

path1=[mfile_path,'/output_s1/pre0.mda'];
path2=[mfile_path,'/output_r1/pre0.mda'];
A1=readmda(path1);
A2=readmda(path2);
A=cat(2,A1,A2);

path0=[mfile_path,'/output_s1r1'];
if (~exist(path0,'dir')) mkdir(path0); end;
writemda(A,[path0,'/pre0.mda']);

end