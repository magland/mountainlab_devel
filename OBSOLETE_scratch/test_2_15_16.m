function test_2_15_16

close all;

path0='franklab/test_hippocampal_02_10_2016/output_tetrode1';
clips=readmda([path0,'/clips.mda']);
clusters=readmda([path0,'/clusters.mda']);
templates=readmda([path0,'/templates.mda']);

view_params.raw=[path0,'/pre2.mda'];
view_params.clusters=[path0,'/clusters.mda'];
view_params.mode='overview2';
ms_mountainview(view_params);

return;

path0='franklab/test_hippocampal_02_10_2016/output_tetrode1';
clips=readmda([path0,'/clips.mda']);
clusters=readmda([path0,'/clusters_split.mda']);
templates=readmda([path0,'/templates_split.mda']);

times=clusters(2,:);
labels=clusters(3,:);
K=max(labels);

k=2;
inds_k=find(labels==k);
clips_k=clips(:,:,inds_k);
clips_k_rand=clips_k(:,:,randsample(length(inds_k),30));

template=templates(:,:,k);
figure; ms_view_templates(template);
figure; ms_view_templates(clips_k_rand);

weights=get_template_weights(template,5);
figure; ms_view_templates(weights);

clips=clips_k_rand;
[M,T,L]=size(clips);
Vclips=reshape(clips,M*T,L);
Vtemplate=reshape(template,M*T,1);
Vweights=reshape(weights,M*T,1);
Vclips=Vclips.*repmat(Vweights,1,L);
Vtemplate=Vtemplate.*Vweights;
Vresid=Vclips-repmat(Vtemplate,1,L);
Vresid_norm=sqrt(sum(Vresid.^2,1));
%Vtemplate_norm=sqrt(Vtemplate'*Vtemplate);
Vweights_norm=sqrt(Vweights'*Vweights);

sigma=1;
expected_vresid_norm=Vweights_norm*sigma;
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

function Y=get_template_weights(template,num_pix)
[M,T]=size(template);

aa=ifftshift(-floor(T/2):-floor(T/2)+T-1);
sig=num_pix;
kernel=exp(-0.5*aa.^2/sig^2);

fhat=fft(abs(template),[],2);
fhat=fhat.*repmat(kernel,M,1);
Y=real(ifft(fhat,[],2));
end

function stdevs=compute_cluster_stats(clips)
stdevs=sqrt(var(clips,3));
end

