function mscmd_isobranch(input_clips_path,output_labels_path,opts)

if (nargin<3) opts=struct; end;
if (~isfield(opts,'isocut_threshold')) opts.isocut_threshold=1.2; end;
if (~isfield(opts,'K_init')) opts.K_init=30; end;
if (~isfield(opts,'num_features')) opts.num_features=3; end;
if (~isfield(opts,'min_cluster_size')) opts.min_cluster_size=500; end;

cmd=sprintf('%s isobranch --input_clips=%s --output_labels=%s --isocut_threshold=%g --K_init=%d --min_cluster_size=%d --num_features=%d ',mscmd_exe,...
    input_clips_path,output_labels_path,opts.isocut_threshold,opts.K_init,opts.min_cluster_size,opts.num_features);

fprintf('\n*** ISOBRANCH ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end