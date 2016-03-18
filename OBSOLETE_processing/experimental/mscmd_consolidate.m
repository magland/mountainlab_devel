function mscmd_consolidate(cluster_path,templates_path,cluster_out_path,templates_out_path,load_channels_path,opts)
% This function will need to be revisited

if (nargin<6) opts=struct; end;

cmd=sprintf('%s consolidate --clusters=%s --templates=%s --cluster_out=%s --templates_out=%s --load_channels_out=%s --coincidence_threshold=%g ',mscmd_exe,cluster_path,templates_path,cluster_out_path,templates_out_path,load_channels_path,opts.coincidence_threshold);

fprintf('\n*** CONSOLIDATE ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end