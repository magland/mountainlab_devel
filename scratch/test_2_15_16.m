function test_2_15_16

path0='franklab/test_hippocampal_02_10_2016/output_tetrode1';
clips=readmda([path0,'/clips.mda']);
clusters=readmda([path0,'/clusters.mda']);

times=clusters(2,:);
labels=clusters(3,:);
K=max(labels);

% stdevs=zeros(size(clips));
% for k=1:K
%     inds_k=find(labels==k);
%     clips_k=clips(:,:,k);
%     stdevs(:,:,k)=compute_cluster_stats(clips_k);
% end;

figure; 
ms_view_templates_from_clips(clips,labels,struct('show_stdev',1));

end

function stdevs=compute_cluster_stats(clips)
stdevs=sqrt(var(clips,3));
end

