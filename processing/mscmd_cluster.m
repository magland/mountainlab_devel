function mscmd_cluster(input_path,output_path,opts)

if (nargin<3) opts=struct; end;

cmd=sprintf('%s cluster --input=%s --output=%s ',mscmd_exe,input_path,output_path);

fprintf('\n*** CLUSTER ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end