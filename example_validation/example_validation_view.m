function example_validation_view

mfile_path=fileparts(mfilename('fullpath'));
path0=sprintf('%s/../example_data',mfile_path);
path1=sprintf('%s/../example_data_reversed',mfile_path);
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s --mode=compare_labels ',exe_fname)];

cmd=[cmd,sprintf('--raw=%s/filt2_white.mda ',path0)];
cmd=[cmd,sprintf('--cluster=%s/cluster0b.mda ',path0)];
cmd=[cmd,sprintf('--cluster2=%s/cluster0b_mapped.mda ',path1)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));

C1=readmda(sprintf('%s/cluster0b.mda',path0));
C2=readmda(sprintf('%s/cluster0b_mapped.mda',path1));

CM=confusion_matrix(C1(2,:),C1(3,:),C2(2,:),C2(3,:));
figure; imagesc(CM');

end

