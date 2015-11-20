function mscmd_consolidate(cluster_path,templates_path,cluster_out_path,templates_out_path,load_channels_path,opts)

if (nargin<5) opts=struct; end;

cmd=sprintf('%s consolidate --cluster=%s --templates=%s --cluster_out=%s --templates_out=%s --load_channels_out=%s ',mscmd_exe,cluster_path,templates_path,cluster_out_path,templates_out_path,load_channels_path);

fprintf('\n*** CONSOLIDATE ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end