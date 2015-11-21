function mscmd_split_clusters(input_path,cluster_path,output_path,opts)

if (nargin<3) opts=struct; end;

cmd=sprintf('%s split_clusters --input=%s --cluster=%s --output=%s --num_features=%d --clip_size=%d ',mscmd_exe,input_path,cluster_path,output_path,opts.num_features,opts.clip_size);

fprintf('\n*** SPLIT CLUSTERS ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end