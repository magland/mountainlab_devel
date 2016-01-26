function mscmd_split_clusters(input_path,cluster_path,output_path,opts)

if (nargin<3) opts=struct; end;
if (~isfield(opts,'ks_threshold')) opts.ks_threshold=1.2; end;
if (~isfield(opts,'K_init')) opts.K_init=30; end;

cmd=sprintf('%s split_clusters --input=%s --cluster=%s --output=%s --num_features=%d --clip_size=%d --ks_threshold=%g ',mscmd_exe,input_path,cluster_path,output_path,opts.num_features,opts.clip_size,opts.ks_threshold);

fprintf('\n*** SPLIT CLUSTERS ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end