function ms_mountainview(view_params)

mfile_path=fileparts(mfilename('fullpath'));
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s ',exe_fname)];

if isfield(view_params,'raw')
    cmd=[cmd,sprintf('--raw=%s ',view_params.raw)];
end;
if isfield(view_params,'clusters')
    cmd=[cmd,sprintf('--clusters=%s ',view_params.clusters)];
end;
if isfield(view_params,'templates')
    cmd=[cmd,sprintf('--templates=%s ',view_params.templates)];
end;
if (isfield(view_params,'clips'))&&(isfield(view_params,'clips_index'))
    cmd=[cmd,sprintf('--clips=%s ',view_params.clips)];
    cmd=[cmd,sprintf('--clips_index=%s ',view_params.clips_index)];
end;
if isfield(view_params,'cross_correlograms')
    cmd=[cmd,sprintf('--cross_correlograms=%s ',view_params.cross_correlograms)];
end;
if isfield(view_params,'locations')
    cmd=[cmd,sprintf('--locations=%s ',view_params.locations)];
end;

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));

end
