function mscmd_isobranch(input_clips_path,output_labels_path,opts)

if (nargin<3) opts=struct; end;
if (~isfield(opts,'isocut_threshold')) opts.isocut_threshold=1.2; end;
if (~isfield(opts,'K_init')) opts.K_init=30; end;
if (~isfield(opts,'num_features')) opts.num_features=6; end;
if (~isfield(opts,'branch_thresholds')) error('Missing required parameter: branch_thresholds'); end;

branch_thresholds_str=array2str(opts.branch_thresholds);

cmd=sprintf('%s isobranch --input_clips=%s --output_labels=%s --branch_thresholds=%s --isocut_threshold=%g --K_init=%d --num_features=%d ',mscmd_exe,...
    input_clips_path,output_labels_path,branch_thresholds_str,opts.isocut_threshold,opts.K_init,opts.num_features);

fprintf('\n*** ISOBRANCH ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end

function str=array2str(X)
str='';
for j=1:length(X)
    str=[str,sprintf('%g',X(j))];
    if (j+1<=length(X)) str=[str,',']; end;
end;
end