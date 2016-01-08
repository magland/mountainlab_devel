function mountainsort(input_file_path,output_dir_path,opts)

total_timer=tic;

basepath=fileparts(mfilename('fullpath'));
addpath([basepath,'/processing']);
addpath([basepath,'/isosplit']);
addpath([basepath,'/msutils']);

locations=opts.locations;
adjacency_matrix=ms_adjacency_matrix(locations,opts.adjacency_radius);

if isfield(opts,'preprocessed_input_file_path')
    fprintf('Reading preprocessed file... '); timer0=tic;
    X=readmda_data_beginning(opts.preprocessed_input_file_path,length(opts.timepoints));
    fprintf('%.2f seconds...\n',toc(timer0));  timer0=tic;
else
    fprintf('Reading... '); timer0=tic;
    X=ms_readdat(input_file_path,opts);
    fprintf('%.2f seconds...\n',toc(timer0));  timer0=tic;
    fprintf('Filtering... ');
    X=ms_filter(X,opts);
    fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;
    fprintf('Normalizing... ');
    X=ms_normalize(X);
    fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;
end;
M=size(X,1);
X0=X;
if (opts.whiten)
%     fprintf('Prewhitening...\n');
%     X=ms_prewhiten(X,opts);
    fprintf('Whitening... ');
    X=ms_whiten(X,opts);
    fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;
    fprintf('Normalizing... ');
    X=ms_normalize(X);
    fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;
end;

opts2=opts; opts2.freq_min=600; opts2.freq_max=2000;
fprintf('Filtering... ');
Y=ms_filter(X,opts2);
fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;
fprintf('Normalizing... ');
Y=ms_normalize(Y);
fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;
all_times=[]; all_labels=[]; primary_channels=[];
current_label=1;
for ch=1:M
    fprintf('Channel %d/%d... ',ch,M);
    [times_pos,times_neg]=ms_detect(Y(ch,:),opts);
    AM=adjacency_matrix; AM(ch,ch)=0;
    neighbor_channels=[ch,find(AM(ch,:))];
    clips_pos=ms_extract_clips(X(neighbor_channels,:),times_pos,opts.clip_size);
    clips_neg=ms_extract_clips(X(neighbor_channels,:),times_neg,opts.clip_size);
    FF_pos=ms_event_features(clips_pos,opts.num_pca_features);
    FF_neg=ms_event_features(clips_neg,opts.num_pca_features);
%     ss_view_clusters(FF_pos,ones(1,size(FF_pos,2))); title(sprintf('%d pos',ch));
%     ss_view_clusters(FF_neg,ones(1,size(FF_neg,2))); title(sprintf('%d neg',ch));
    labels_pos=ms_cluster(FF_pos); Kpos=max(labels_pos); if (isempty(Kpos)) Kpos=0; end;
    labels_neg=ms_cluster(FF_neg); Kneg=max(labels_neg); if (isempty(Kneg)) Kneg=0; end;
    times0=[times_pos,times_neg];
    labels0=[labels_pos,labels_neg+Kpos];
    clips0=ms_extract_clips(X,times0,opts.clip_size); %Use X or X0 here?
    fprintf('spl ');
    labels0=split_clusters(clips0,labels0,3);
    fprintf('rem ');
    [times0,labels0]=remove_outliers(times0,labels0,clips0,opts.num_pca_features,opts.cluster_outlier_alpha);
    templates0=compute_templates(X,times0,labels0,opts.clip_size); %Use X or X0 here?
    fprintf('con ');
    K_before_consolidate=max(labels0);
    [times0,labels0,templates0]=consolidate_clusters(times0,labels0,templates0,ch,opts.min_cluster_size);
    
    K0=max(labels0);
    if (K0>0)
        labels0=labels0+current_label-1;
        current_label=current_label+K0;
        all_times=[all_times,times0];
        all_labels=[all_labels,labels0];
        primary_channels=[primary_channels,ones(1,K0)*ch];
    end;
    
    fprintf('%d clusters, using %d (%d events)... ',K_before_consolidate,K0,length(times0));
    fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;
end;
fprintf('Sorting times/labels... ');
[all_times,sort_inds]=sort(all_times);
all_labels=all_labels(sort_inds);
fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;

fprintf('***** Detected %d events and %d neurons... *****\n',length(all_times),max(all_labels));

fprintf('Computing cross correlograms...\n');
[CC,CCmda]=ms_cross_correlograms(all_times,all_labels,opts.cross_correlogram_max_dt);
fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;

fprintf('writing output... ');
outdir=output_dir_path;
if (~exist(outdir,'dir'))
    mkdir(outdir);
end;

if (opts.whiten)
    writemda(X,[outdir,'/raw_white.mda']);
    writemda(compute_templates(X,all_times,all_labels,opts.clip_size),[outdir,'/templates_white.mda']);
else
    if (exist([outdir,'/raw_white.mda'],'file'))
        delete([outdir,'/raw_white.mda']);
    end;
    if (exist([outdir,'/templates_white.mda'],'file'))
        delete([outdir,'/templates_white.mda']);
    end;
end
writemda(X0,[outdir,'/raw.mda']);
writemda(all_times,[outdir,'/times.mda']);
writemda(all_labels,[outdir,'/labels.mda']);
writemda(compute_templates(X0,all_times,all_labels,opts.clip_size),[outdir,'/templates.mda']);
writemda(adjacency_matrix,[outdir,'/adjacency.mda']);
writemda(locations,[outdir,'/locations.mda']);
writemda(primary_channels,[outdir,'/primary_channels.mda']);
writemda(CCmda,[outdir,'/cross-correlograms.mda']);
fprintf('%.2f seconds...\n',toc(timer0)); timer0=tic;

fprintf('Total time for mountainsort: %.2f\n',toc(total_timer));

end

function labels2=split_clusters(X,labels,num_pca_features)
K=max(labels);
labels2=zeros(size(labels));
current_label=1;
for k=1:K
    fprintf(':');
    inds=find(labels==k);
    X0=X(:,:,inds);
    FF0=ms_event_features(X0,num_pca_features);
    labels000=ms_cluster(FF0);
    K000=max(labels000);
    if (K000>1)
        fprintf('(%d)',K000);
    end;
    for jj=1:K000
        indsjj=find(labels000==jj);
        labels2(inds(indsjj))=current_label;
        current_label=current_label+1;
    end;
end;
fprintf(' ');
end

function [times,labels,templates]=consolidate_clusters(times,labels,templates,k,min_cluster_size)
sizes=squeeze(sum(templates.^2,2));
max_sizes=(max(sizes,[],1));
rel_sizes=sizes(k,:)./max_sizes;
to_use=(rel_sizes>=0.9);
K=max(labels);

for k=1:K
    ct=length(find(labels==k));
    if (ct<min_cluster_size)
        to_use(k)=0;
    end;
end;

labels_to_use=find(to_use);
ret=length(labels_to_use);

templates=templates(:,:,labels_to_use);

times2=[];
labels2=[];
for ii=1:length(labels_to_use)

    inds0=find(labels==labels_to_use(ii));
    times2=[times2,times(inds0)];
    labels2=[labels2,ones(size(times(inds0)))*ii];
end;
[times,sort_inds]=sort(times2);
labels=labels2(sort_inds);

end

function templates=compute_templates(X,times,labels,clip_size)
[M,N]=size(X);
T=clip_size;
clips=ms_extract_clips(X,times,clip_size);
K=max(labels);
templates=zeros(M,T,K);
for k=1:K
    inds=find(labels==k);
    templates(:,:,k)=median(clips(:,:,inds),3);
end;
end

function [times2,labels2]=remove_outliers(times,labels,X,num_pca_features,cluster_outlier_alpha)

K=max(labels);
use_it=ones(size(labels));
for k=1:K
    inds=find(labels==k);
    FF=ms_event_features(X(:,:,inds),num_pca_features);
    inds2=find_outliers(FF,cluster_outlier_alpha);
    use_it(inds(inds2))=0;
end;

ok_inds=find(use_it);
times2=times(ok_inds);
labels2=labels(ok_inds);

end

function inds=find_outliers(X,alpha)

X0=median(X,2);
X=X-repmat(X0,1,size(X,2));
for j=1:size(X,1)
    X(j,:)=X(j,:)/sqrt(var(X(j,:)));
end;

vals=sum(X.^2,1);
p=1-chi2cdf(vals,size(X,1));
inds=find(p<alpha);

end

