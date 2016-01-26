function mscmd_cluster(input_path,output_path,opts)

if (nargin<3) opts=struct; end;
if (~isfield(opts,'ks_threshold')) opts.ks_threshold=1.2; end;
if (~isfield(opts,'K_init')) opts.K_init=30; end;

cmd=sprintf('%s cluster --input=%s --output=%s --ks_threshold=%g --K_init=%d ',mscmd_exe,input_path,output_path,opts.ks_threshold,opts.K_init);

fprintf('\n*** CLUSTER ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end