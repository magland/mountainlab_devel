function test_2016_01_21_validate_view

mfile_path=fileparts(mfilename('fullpath'));
addpath(sprintf('%s/..',mfile_path));
outputdir_path=sprintf('%s/output',mfile_path);
path0=sprintf('%s/../experiment1/output',mfile_path);

raw=sprintf('%s/pre1.mda',path0);
clusters1=sprintf('%s/clusters.mda',path0);
clusters2=sprintf('%s/clusters.mda',outputdir_path);

open_mv_compare(raw,clusters1,clusters2);

VM=readmda(sprintf('%s/validation_matrix.mda',outputdir_path));
[K1,K2]=size(VM); K1=K1-1; K2=K2-1;
VM=normalize_confusion_matrix(VM);
reordering=get_best_reordering(VM');
VM=VM(:,[reordering,K2+1]);
figure; imagesc(VM'); colorbar; xlabel('Original Run'); ylabel('Validation Run');
title('Validation Matrix');

figure;
templates=readmda(sprintf('%s/templates_raw.mda',path0));
ms_view_templates(templates);

end

function CM=normalize_confusion_matrix(CM)
[K1,K2]=size(CM); K1=K1-1; K2=K2-1;
CM_row_normalized=CM./repmat(sum(CM,2),1,K2+1);
CM_column_normalized=CM./repmat(sum(CM,1),K1+1,1);
CM_normalized=(CM_row_normalized+CM_column_normalized)/2;
CM=CM_normalized;
end

function open_mv_compare(raw,clusters1,clusters2)
mfile_path=fileparts(mfilename('fullpath'));
exe_fname=sprintf('%s/../../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s --mode=compare_labels ',exe_fname)];

cmd=[cmd,sprintf('--raw=%s ',raw)];
cmd=[cmd,sprintf('--cluster=%s ',clusters1)];
cmd=[cmd,sprintf('--cluster2=%s ',clusters2)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));
end

function reordering=get_best_reordering(CM)
[K1,K2]=size(CM);
K1=K1-1;
K2=K2-1;
HH=Hungarian(-CM(1:K1,1:K2));

map12=zeros(1,K1);
for j=1:K1
    ind0=find(HH(j,1:K2)==1);
    if (length(ind0)==0) map12(j)=0;
    else map12(j)=ind0(1);
    end;
end;

reordering=map12;

used=zeros(1,K1);
used(reordering(reordering>0))=1;
unused_inds=find(used==0);
inds=find(reordering==0);
reordering(inds)=unused_inds;

[~,reordering]=sort(reordering);

end
