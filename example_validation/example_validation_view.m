function example_validation_view

mfile_path=fileparts(mfilename('fullpath'));
path0=sprintf('%s/output',mfile_path);
path1=sprintf('%s/output_reversed',mfile_path);
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s --mode=compare_labels ',exe_fname)];

cmd=[cmd,sprintf('--raw=%s/filt2_white.mda ',path0)];
cmd=[cmd,sprintf('--cluster=%s/clusters.mda ',path0)];
cmd=[cmd,sprintf('--cluster2=%s/clusters_mapped.mda ',path1)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));

C1=readmda(sprintf('%s/clusters.mda',path0));
C2=readmda(sprintf('%s/clusters_mapped.mda',path1));

CM=confusion_matrix(C1(2,:),C1(3,:),C2(2,:),C2(3,:));
figure; imagesc(CM');

[K1,K2]=size(CM);

CM2=CM;
for k1=1:K1
    for k2=1:K2
        CM2(k1,k2)=2*CM(k1,k2)/(sum(CM(k1,:))+sum(CM(:,k2)));
    end;
end;
figure; imagesc(CM2');

view_output(path0,1);

end

function view_output(output_path,whitened)

mfile_path=fileparts(mfilename('fullpath'));
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s ',exe_fname)];
if (whitened)
    cmd=[cmd,sprintf('--raw=%s/filt2_white.mda ',output_path)];
    cmd=[cmd,sprintf('--templates=%s/templates0_filt2_white.mda ',output_path)];
    cmd=[cmd,sprintf('--clips=%s/clips_filt2_white.mda ',output_path)];
else
    cmd=[cmd,sprintf('--raw=%s/filt.mda ',output_path)];
    cmd=[cmd,sprintf('--templates=%s/templates0_raw.mda ',output_path)];
    cmd=[cmd,sprintf('--clips=%s/clips_filt.mda ',output_path)];
end;

cmd=[cmd,sprintf('--primary-channels=%s/load_channels0.mda ',output_path)];

cmd=[cmd,sprintf('--cluster=%s/clusters.mda ',output_path)];
cmd=[cmd,sprintf('--locations=%s/locations.mda ',output_path)];
cmd=[cmd,sprintf('--cross-correlograms=%s/cross-correlograms.mda ',output_path)];

cmd=[cmd,sprintf('--clips-index=%s/clips_filt_index.mda',output_path)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));

end

