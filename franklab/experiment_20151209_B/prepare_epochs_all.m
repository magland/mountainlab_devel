function prepare_epochs_all

mfile_path=fileparts(mfilename('fullpath'));

strings={'s1','r1','s2','r2','s3','r3'};
for j=1:length(strings)
    str=strings{j};
    disp(str);
    path=[mfile_path,sprintf('/output_%s/pre0.mda',str)];
    A=readmda(path);
    if (j==1) AA=A;
    else AA=cat(2,AA,zeros(size(A,1),10e6),A);
    end;
end;

path0=[mfile_path,'/output_all_epochs'];
if (~exist(path0,'dir')) mkdir(path0); end;
writemda32(AA,[path0,'/pre0.mda']);

end