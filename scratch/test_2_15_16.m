function test_2_15_16

path0='franklab/test_hippocampal_02_10_2016/output_tetrode1';
clips=readmda([path0,'/clips.mda']);
clusters=readmda([path0,'/clusters_split.mda']);
templates=readmda([path0,'/templates_split.mda']);

times=clusters(2,:);
labels=clusters(3,:);
K=max(labels);

k=1;
inds_k=find(labels==k);
clips_k=clips(:,:,inds_k);
clips_k_rand=clips_k(:,:,randsample(length(inds_k),30));
figure; ms_view_templates(cat(3,templates(:,:,k),clips_k_rand));

template=templates(:,:,k);
clips=clips_k_rand;
[M,T,L]=size(clips);
Vclips=reshape(clips,M*T,L);
Vtemplate=reshape(template,M*T,1);
Vclips=Vclips.*abs(repmat(Vtemplate,1,L));
Vtemplate=Vtemplate.*abs(Vtemplate);
Vresid=Vclips-repmat(Vtemplate,1,L);
Vresid_norm=sqrt(sum(Vresid.^2,1));
Vtemplate_norm=sqrt(Vtemplate'*Vtemplate);

sigma=1;
expected_vresid_norm=Vtemplate_norm*sigma;
score=Vresid_norm./expected_vresid_norm;

disp(score);

% stdevs=zeros(size(clips));
% for k=1:K
%     inds_k=find(labels==k);
%     clips_k=clips(:,:,k);
%     stdevs(:,:,k)=compute_cluster_stats(clips_k);
% end;

%figure; 
%ms_view_templates_from_clips(clips,labels,struct('show_stdev',1));

end

function stdevs=compute_cluster_stats(clips)
stdevs=sqrt(var(clips,3));
end

