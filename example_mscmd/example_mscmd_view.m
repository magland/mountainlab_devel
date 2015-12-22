function example_mscmd_view(whitened)

if nargin<1
    whitened=input('Enter 1 for whitened 0 for non-whitened view: ');
end

mfile_path=fileparts(mfilename('fullpath'));
path=sprintf('%s/../example_data',mfile_path);
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s ',exe_fname)];
if whitened
    cmd=[cmd,sprintf('--raw=%s/filt2_white.mda ',path)];
    cmd=[cmd,sprintf('--templates=%s/templates0_filt2_white.mda ',path)];
    cmd=[cmd,sprintf('--clips=%s/clips_filt2_white.mda ',path)];
else
    cmd=[cmd,sprintf('--raw=%s/filt.mda ',path)];
    cmd=[cmd,sprintf('--templates=%s/templates0_raw.mda ',path)];
    cmd=[cmd,sprintf('--clips=%s/clips_filt.mda ',path)];
end;
cmd=[cmd,sprintf('--cluster=%s/cluster0b.mda --primary-channels=%s/load_channels0.mda ',path,path)];
cmd=[cmd,sprintf('--locations=%s/locations.mda --clips-index=%s/clips_filt_index.mda --cross-correlograms=%s/cross-correlograms.mda ',path,path,path)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));

end

