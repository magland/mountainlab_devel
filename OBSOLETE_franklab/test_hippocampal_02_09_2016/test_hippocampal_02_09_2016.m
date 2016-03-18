function test_hippocampal_02_09_2016

close all; drawnow;

%%%% Parameters and settings
tetrode_num=1;
plausibility_threshold=0.6;
merge_threshold=0.9; %Keep this high for now
tt_range=[5,15];
num_tt_steps=10;
tt_overlap=1;
num_features=6;
cross_correlograms_max_dt=6000;
sigma=1.5;
o_mild_filter.samplefreq=30000;
o_mild_filter.freq_min=50;
o_mild_filter.freq_max=10000;
o_mild_filter.outlier_threshold=500;
o_filter.samplefreq=30000;
o_filter.freq_min=100;
o_filter.freq_max=10000;
o_filter.outlier_threshold=500;
o_detect.threshold=tt_range(1);
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=15;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=120;
o_whiten=struct;

%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode'];
path0=[mfile_path,sprintf('/output_tetrode%d',tetrode_num)];
if ~exist(path0,'dir') mkdir(path0); end;

%%%% Extract raw data
extract_raw_data(raw_path,path0,tetrode_num);

%%%% Preprocessing
mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre1.mda'],o_filter);
mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre0_mild.mda'],o_mild_filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],o_whiten);
mscmd_detect([path0,'/pre2.mda'],[path0,'/detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[path0,'/detect.mda'],[path0,'/clips.mda'],o_extract_clips);

%%%% Reading clips
fprintf('Reading...\n');
clips=readmda([path0,'/clips.mda']);
[M,T,NC]=size(clips);

%%%% Shell cluster
fprintf('Shell cluster...\n');
[templates0,clusterings,peaks0]=shell_cluster(clips,tt_range,num_tt_steps,tt_overlap,num_features);
all_templates=group_templates(templates0);
figure; ms_view_templates(all_templates);
title('Unmerged templates'); drawnow;
[~,~,K_all]=size(all_templates);
title('All templates'); drawnow;

%%%% Posterior probabilities and plausibility factors
fprintf('Posterior probabilities and plausibility factors...');
[probs,plausibility_factors]=compute_posterior_probabilities(clips,all_templates,sigma);
corr_matrix=corrcoef(probs');
figure; imagesc(corr_matrix');
title('Correlation matrix for posterior classification prob vectors'); drawnow;
figure; hist(max(plausibility_factors,[],1),linspace(-0.02,1.02,1000));
title('Plausibility factors for all events'); drawnow;

%%%% Times and labels
fprintf('Computing times/labels and coincidence matrix...\n');
detect=readmda([path0,'/detect.mda']);
times0=detect(2,:);
times=[]; labels_unmerged=[]; peaks=[]; used=zeros(size(times0));
[~,best_ks]=max(probs,[],1);
for k=1:size(probs,1);
    inds_k=find((plausibility_factors(k,:)>=plausibility_threshold)&(best_ks==k));
    times=[times,times0(inds_k)];
    labels_unmerged=[labels_unmerged,ones(1,length(inds_k))*k];
    peaks=[peaks,peaks0(inds_k)];
    used(inds_k)=1;
end;
inds_0=find(used==0);
times=[times,times0(inds_0)];
labels_unmerged=[labels_unmerged,ones(1,length(inds_0))*0];
peaks=[peaks,peaks0(inds_0)];
[~,sort_inds]=sort(times); times=times(sort_inds); labels_unmerged=labels_unmerged(sort_inds); peaks=peaks(sort_inds);
clusters_unmerged=zeros(4,length(times));
clusters_unmerged(2,:)=times;
clusters_unmerged(3,:)=labels_unmerged;
clusters_unmerged(4,:)=peaks;

%%%% Merge
cliques=ms_greedy_cliques(corr_matrix>=merge_threshold);
label_map=ones(1,max(labels_unmerged));
for jj=1:length(cliques)
    label_map(cliques{jj})=jj;
end;
clusters=clusters_unmerged;
nonzero_inds=find(clusters_unmerged(3,:)>0);
clusters(3,nonzero_inds)=label_map(clusters(3,nonzero_inds));
labels=clusters(3,:);

%%%% Writing output and preparing view
fprintf('Writing output and preparing view...\n');
pre2=readmda([path0,'/pre2.mda']);
[clips1,clips1_index]=ms_create_clips_index(ms_extract_clips(pre2,times,o_extract_clips.clip_size),labels);
writemda(clips1,[path0,'/clips0.mda']);
writemda(clips1_index,[path0,'/clips0_index.mda']);
writemda(clusters,[path0,'/clusters.mda']);
writemda(corr_matrix,[path0,'/correlation_matrix.mda']);

%%%% Cross correlograms and templates
mscmd_cross_correlograms([path0,'/clusters.mda'],[path0,'/cross_correlograms.mda'],cross_correlograms_max_dt);
mscmd_templates([path0,'/pre0_mild.mda'],[path0,'/clusters.mda'],[path0,'/templates_raw.mda'],struct('clip_size',200));
mscmd_templates([path0,'/pre2.mda'],[path0,'/clusters.mda'],[path0,'/templates.mda'],struct('clip_size',200));
templates=readmda([path0,'/templates.mda']);
figure; ms_view_templates(templates);
title('Merged templates'); drawnow;

%%%% MountainView
view_params.raw=[path0,'/pre2.mda'];
view_params.clusters=[path0,'/clusters.mda'];
view_params.cross_correlograms=[path0,'/cross_correlograms.mda'];
view_params.templates=[path0,'/templates.mda'];
view_params.clips=[path0,'/clips0.mda'];
view_params.clips_index=[path0,'/clips0_index.mda'];
ms_mountainview(view_params);

end

function [probs,plausibility_factors]=compute_posterior_probabilities(clips,templates,sigma)
fprintf('Computing ips... ');
[M,T,K]=size(templates);
[~,~,NC]=size(clips);
ips=zeros(K,NC);
for k=1:K
    fprintf('%d ',k);
    ips(k,:)=squeeze(sum(sum(clips.*repmat(templates(:,:,k),1,1,NC),1),2));
end;
fprintf('\nComputing norms, etc...\n');
template_norms=reshape(sqrt(sum(sum(templates.^2,1),2)),1,K);
clip_norms=reshape(sqrt(sum(sum(clips.^2,1),2)),1,NC);
diffsqr=repmat(clip_norms,K,1).^2-2*ips+repmat(template_norms',1,NC).^2;

%logprobs1=-2*diffsqr/sigma^2/(M*T);
logprobs1=log(get_pvals(sqrt(diffsqr),sqrt(M*T)));
logmargins=logsumexp(cat(1,logprobs1,log(0.2*ones(1,NC))),1);
%logmargins=logsumexp(logprobs1);
logprobs=logprobs1-repmat(logmargins,K,1);
probs=exp(logprobs);

if nargout>1
    % K x K
    normalization_probs=compute_posterior_probabilities(templates,templates,sigma);
    % K x 1
    normalization_probs=diag(normalization_probs);
    % K x NC
    plausibility_factors=min(1,probs./repmat(normalization_probs,1,NC));
end;

end

function pvals=get_pvals(X,mu)
pvals=min(1,(mu./X).^2);
end



function [templates,clusterings,clip_peaks]=shell_cluster(clips,tt_range,num_tt_steps,tt_overlap,num_features)
[M,T,NC]=size(clips);

fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

sorted_peaks=sort(clip_peaks);
inds00=find((sorted_peaks>=tt_range(1))&(sorted_peaks<=tt_range(2)));
sorted_peaks=sorted_peaks(inds00);
tt1_list=zeros(1,num_tt_steps);
tt2_list=zeros(1,num_tt_steps);
for jj=1:num_tt_steps
    ind1=max(1,ceil(length(sorted_peaks)*(jj-1)/num_tt_steps+1));
    ind2=min(length(sorted_peaks),ceil(length(sorted_peaks)*(jj-1+1+tt_overlap)/num_tt_steps+1));
    tt1_list(jj)=sorted_peaks(ind1);
    tt2_list(jj)=sorted_peaks(ind2);
end;

clusterings={};
for ii=1:length(tt1_list)
    tt1=tt1_list(ii);
    tt2=tt2_list(ii);
    if (ii==length(tt1_list)), tt2=inf; end;
    fprintf('tt1=%g tt2=%g... ',tt1,tt2);
    inds_tt=find((clip_peaks>=tt1)&(clip_peaks<=tt2));
    CC.tt=tt1;
    CC.inds=inds_tt;
    if (length(inds_tt)>1)
        clips_tt=clips(:,:,inds_tt);
        fprintf('features... ');
        [FF_tt,subspace_tt]=ms_event_features(clips_tt,num_features);
        fprintf('isosplit... ');
        labels_tt=isosplit(FF_tt);
        K=max(labels_tt);
        fprintf('K=%d\n',K);
        CC.labels=labels_tt;
        CC.K=K;
        CC.features=FF_tt;
        CC.clips=clips_tt;
        CC.medoids=zeros(M,T,K);
        for k=1:K
            FF_tt_k=FF_tt(:,find(labels_tt==k));
            FF_tt_k_medoid=compute_medoid(FF_tt_k);
            template0=zeros(M,T);
            for jjj=1:size(subspace_tt,3)
                template0=template0+FF_tt_k_medoid(jjj)*subspace_tt(:,:,jjj);
            end;
            CC.medoids(:,:,k)=template0;
        end;
    else
        CC.labels=[];
        CC.K=0;
        CC.features=zeros(num_features,0);
        CC.clips=zeros(M,T,0);
        CC.medoids=zeros(M,T,0);
    end;
    clusterings{end+1}=CC;
end

for ii=1:length(clusterings)-1
    K1=clusterings{ii}.K;
    K2=clusterings{ii+1}.K;
    LM=zeros(K1,K2);
    [inds_intersect,ii1,ii2]=intersect(clusterings{ii}.inds,clusterings{ii+1}.inds);
    if (length(inds_intersect)>0)
        labels1=clusterings{ii}.labels(ii1);
        labels2=clusterings{ii+1}.labels(ii2);
        for k1=1:K1
            for k2=1:K2
                LM(k1,k2)=length(find((labels1==k1)&(labels2==k2)));
            end;
        end;
    end;
    clusterings{ii}.label_matches=LM;
end;
clusterings{end}.label_matches=zeros(clusterings{end}.K,0);

templates=zeros(M,T,0);
for ii=1:length(clusterings)
    templates=cat(3,templates,clusterings{ii}.medoids);
end;

end

function extract_raw_data(raw_path,output_path,tetrode_num)

raw_mat_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mat',raw_path);
raw_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mda',raw_path);
tetrode_fname=sprintf('%s/pre0.mda',output_path);

if (~exist(raw_mda_fname,'file'))
    fprintf('Loading raw data...\n');
    L=load(raw_mat_fname);
    raw=L.dl12_20151208_NNF_r1_tet16_17.channelData';
    fprintf('Writing raw data...\n');
    writemda(raw,raw_mda_fname);
end;

if (~exist(tetrode_fname,'file'))
    fprintf('Reading raw data...\n');
    raw=readmda(raw_mda_fname);

    if (tetrode_num==1)
        tetrode=raw([1,3:6],(1e6+1):26e6);
        tetrode=tetrode(2:end,:)-repmat(tetrode(1,:),size(tetrode,1)-1,1);
    elseif (tetrode_num==2)
        tetrode=raw([1,7:10],(1e6+1):26e6);
        tetrode=tetrode(2:end,:)-repmat(tetrode(1,:),size(tetrode,1)-1,1);
    end;
    fprintf('Writing tetrode data...\n');
    writemda(tetrode,tetrode_fname);
end;

L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];
writemda(L,[output_path,'/locations.mda']);

end

function m = compute_medoid(X)
[M,N]=size(X);
dists=zeros(N,N);
for m=1:M
    [grid1,grid2]=ndgrid(X(m,:),X(m,:));
    dists=dists+sqrt((grid1-grid2).^2);
    %dists=dists+(grid1-grid2).^2;
end;
avg_dists=mean(dists,1);
[~,ind]=min(avg_dists);
m=X(:,ind); 
end

function templates=group_templates(templates)

[M,T,K]=size(templates);
peaks_pos=squeeze(max(templates(:,T/2+1,:),[],1))';
peaks_neg=-squeeze(max(-templates(:,T/2+1,:),[],1))';
peaks=peaks_pos.*(abs(peaks_pos)>abs(peaks_neg))+peaks_neg.*(abs(peaks_pos)<abs(peaks_neg));

[~,sort_inds]=sort(peaks);
templates=templates(:,:,sort_inds);

[M,T,K]=size(templates);
MM=zeros(K,K);
templates=templates-repmat(mean(templates,2),1,T,1);
for k1=1:K
    template1=templates(:,:,k1);
    for k2=1:K
        template2=templates(:,:,k2);
        MM(k1,k2)=compute_correlation(template1(:),template2(:));
    end;
end;

%figure; imagesc(MM'); colorbar;
MM2=MM>0.8;
CC=ms_greedy_cliques(MM2);
%figure; imagesc(MM2'); colorbar;

new_templates=zeros(M,T,0);
for j=1:length(CC)
    new_templates=cat(3,new_templates,templates(:,:,CC{j}));
    %new_templates=cat(3,new_templates,zeros(M,T,2));
end;
templates=new_templates;

end

function ret=compute_correlation(v1,v2)
ret=sum(v1.*v2)/sqrt(sum(v1.^2)*sum(v2.^2));
end
