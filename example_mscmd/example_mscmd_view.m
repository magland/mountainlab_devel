function example_mscmd_view

mfile_path=fileparts(mfilename('fullpath'));
path=sprintf('%s/../example_data',mfile_path);
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s --raw=%s/filt.mda --raw-whitened=%s/filt2_white.mda ',exe_fname,path,path)];
cmd=[cmd,sprintf('--cluster=%s/cluster0.mda --templates=%s/templates0_raw.mda ',path,path)];
cmd=[cmd,sprintf('--templates-whitened=%s/templates0.mda --primary-channels=%s/load_channels0.mda ',path,path)];
cmd=[cmd,sprintf('--locations=%s/locations.mda --cross-correlograms=%s/cross-correlograms.mda ',path,path)];
cmd=[cmd,sprintf('--clips=%s/clips_filt.mda --clips-index=%s/clips_filt_index.mda',path,path)];

fprintf('%s\n',cmd);
system(cmd);

end

