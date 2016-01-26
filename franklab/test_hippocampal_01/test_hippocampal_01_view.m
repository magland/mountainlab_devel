function test_hippocampal_01_view(whitened)

if nargin<1 whitened=0; end;

mfile_path=fileparts(mfilename('fullpath'));
path0=sprintf('%s/output',mfile_path);
exe_fname=sprintf('%s/../../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s ',exe_fname)];

if whitened
    cmd=[cmd,sprintf('--raw=%s/pre3.mda ',path0)];
    cmd=[cmd,sprintf('--templates=%s/templates.mda ',path0)];
    cmd=[cmd,sprintf('--clips=%s/clips.mda ',path0)];
else
    cmd=[cmd,sprintf('--raw=%s/pre1.mda ',path0)];
    cmd=[cmd,sprintf('--templates=%s/templates_pre1.mda ',path0)];
    cmd=[cmd,sprintf('--clips=%s/clips_pre1.mda ',path0)];
end;


cmd=[cmd,sprintf('--cluster=%s/clusters.mda ',path0)];
cmd=[cmd,sprintf('--locations=%s/locations.mda --clips-index=%s/clips_index.mda --cross-correlograms=%s/cross_correlograms.mda ',path0,path0,path0)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));

end
