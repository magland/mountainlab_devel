function step3_cluster(opts)

timerA=tic;

detect_times_prefix=[opts.working_path,'/detect/times_'];
detect_clips_prefix=[opts.working_path,'/detect/clips_'];
cluster_times_prefix=[opts.working_path,'/cluster/times_'];
cluster_labels_prefix=[opts.working_path,'/cluster/labels_'];
cluster_templates_prefix=[opts.working_path,'/cluster/templates_'];
cluster_features_prefix=[opts.working_path,'/cluster/features_'];
adjacency_path=[opts.output_path,'/adjacency.mda'];

AM=readmda(adjacency_path);
M=size(AM,1);

fprintf('Clustering...\n');
for j=1:M
    fname_detect_times_pos=[detect_times_prefix,sprintf('pos_%d.mda',j)];
    fname_detect_times_neg=[detect_times_prefix,sprintf('neg_%d.mda',j)];
    fname_detect_clips_pos=[detect_clips_prefix,sprintf('pos_%d.mda',j)];
    fname_detect_clips_neg=[detect_clips_prefix,sprintf('neg_%d.mda',j)];
    fname_cluster_times=[cluster_times_prefix,sprintf('%d.mda',j)];
    fname_cluster_labels=[cluster_labels_prefix,sprintf('%d.mda',j)];
    fname_cluster_labels_pos=[cluster_labels_prefix,sprintf('pos_%d.mda',j)];
    fname_cluster_labels_neg=[cluster_labels_prefix,sprintf('neg_%d.mda',j)];
    fname_cluster_templates=[cluster_templates_prefix,sprintf('%d.mda',j)];
    fname_cluster_features_pos=[cluster_features_prefix,sprintf('pos_%d.mda',j)];
    fname_cluster_features_neg=[cluster_features_prefix,sprintf('neg_%d.mda',j)];
    if (~exist(fname_cluster_labels,'file'))||(~exist(fname_cluster_templates,'file'))||(~exist(fname_cluster_features_pos,'file'))||(~exist(fname_cluster_features_neg,'file'))
        fprintf('Reading %s and %s...\n',fname_detect_times_pos,fname_detect_times_neg);
        detect_times_pos=readmda(fname_detect_times_pos);
        detect_times_neg=readmda(fname_detect_times_neg);
        fprintf('Reading %s and %s...\n',fname_detect_clips_pos,fname_detect_clips_neg);
        clips_pos=readmda(fname_detect_clips_pos);
        clips_neg=readmda(fname_detect_clips_neg);
        [M,T,Nclips_pos]=size(clips_pos);
        [M,T,Nclips_neg]=size(clips_neg);
        fprintf('Extracting %d features...',opts.num_cluster_features);
        AM(j,j)=0;
        adjacent_channels=[j;find(AM(:,j))];
        %mean_clips_pos=mean(clips_pos(adjacent_channels,:,:),3);
        %mean_clips_neg=mean(clips_neg(adjacent_channels,:,:),3);
        %FF_pos=ss_eventfeatures(clips_pos(adjacent_channels,:,:)-repmat(mean_clips_pos,1,1,Nclips_pos));
        %FF_neg=ss_eventfeatures(clips_neg(adjacent_channels,:,:)-repmat(mean_clips_neg,1,1,Nclips_neg));
        
        %norms_pos=repmat(sqrt(sum(sum(clips_pos(adjacent_channels,:,:).^2,1),2)),length(adjacent_channels),T,1);
        %norms_neg=repmat(sqrt(sum(sum(clips_neg(adjacent_channels,:,:).^2,1),2)),length(adjacent_channels),T,1);
        %clipsFF_pos=clips_pos(adjacent_channels,:,:)./norms_pos;
        %clipsFF_neg=clips_neg(adjacent_channels,:,:)./norms_neg;
        
        clipsFF_pos=clips_pos(adjacent_channels,:,:);
        clipsFF_neg=clips_neg(adjacent_channels,:,:);
        
        FF_pos=ss_eventfeatures(clipsFF_pos);
        FF_neg=ss_eventfeatures(clipsFF_neg);
        FF_pos=FF_pos(1:opts.num_cluster_features,:);
        FF_neg=FF_neg(1:opts.num_cluster_features,:);
        
        fprintf('ISO-SPLIT...\n');
        times=[];
        labels=[];
        if Nclips_pos>0
            labels_pos=do_cluster(clipsFF_pos,opts);
            Kpos=max(labels_pos);
            times=[times,detect_times_pos];
            labels=[labels,labels_pos];
        else
            labels_pos=[];
            Kpos=0;
        end;
        if Nclips_neg>0
            labels_neg=do_cluster(clipsFF_neg,opts);
            times=[times,detect_times_neg];
            labels=[labels,labels_neg+Kpos];
        else
            labels_neg=[];
        end;
        
        clips=cat(3,clips_pos,clips_neg);
        [times,sort_inds]=sort(times);
        labels=labels(sort_inds);
        clips=clips(:,:,sort_inds);
        
        fprintf('Computing average templates... ');
        K=max(labels);
        WF=zeros(M,T,K);
        for k=1:K
            inds=find(labels==k);
            WF(:,:,k)=median(clips(:,:,inds),3);
        end;
        
%         fprintf('Removing the noise clusters...\n');
%         noise_clusters=[];
%         for k=1:K
%             if (is_noise_cluster(
%         end;
        
%         fprintf('Removing the zero clusters...\n');
%         magnitudes=squeeze(sqrt(sum(sum(WF.^2,1),2)));
%         [~,mag_sort_inds]=sort(magnitudes);
%         zero_label_1=mag_sort_inds(1);
%         zero_label_2=mag_sort_inds(2);
%         labels(find(labels==zero_label_1))=0;
%         labels(find(labels==zero_label_2))=0;
%         okay_inds=find(labels>0);
%         times=times(okay_inds);
%         labels=labels(okay_inds);
%         clips=clips(:,:,okay_inds);
%         labels0=labels;
%         labels(find(labels0>zero_label_1))=labels(find(labels0>zero_label_1))-1;
%         labels(find(labels0>zero_label_2))=labels(find(labels0>zero_label_2))-1;
%         if (K-2~=max(labels)) error('Unexpected problem!'); end;
        
        %[times,labels,clips]=remove_small_clusters(times,labels,clips,60);
        
        [times,labels,clips]=refine_clusters(times,labels,clips,adjacent_channels,opts);
        
        K=max(labels);
        
        WF=zeros(M,T,K);
        for k=1:K
            inds=find(labels==k);
            WF(:,:,k)=median(clips(:,:,inds),3);
        end;
        
        K_to_use=compute_number_of_templates_to_use(WF,j);
        
        fprintf('Writing %s... ',fname_cluster_templates);
        writemda(WF,fname_cluster_templates);
        fprintf('Writing %s...\n',fname_cluster_times);
        writemda(times,fname_cluster_times);
        fprintf('Writing %s (%d clusters detected, using %d)...\n',fname_cluster_labels,K,K_to_use);
        writemda(labels,fname_cluster_labels);
        fprintf('Writing %s and %s...\n',fname_cluster_labels_pos,fname_cluster_labels_neg);
        writemda(labels_pos,fname_cluster_labels_pos);
        writemda(labels_neg,fname_cluster_labels_neg);
        fprintf('Writing %s and %s... ',fname_cluster_features_pos,fname_cluster_features_neg);
        writemda(FF_pos,fname_cluster_features_pos);
        writemda(FF_neg,fname_cluster_features_neg);
    else
        fprintf('Files exist %s, %s, %s and %s...\n',fname_cluster_labels,fname_cluster_templates,fname_cluster_features_pos,fname_cluster_features_neg);
    end;
end;

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end

function [times2,labels2,clips2]=remove_small_clusters(times,labels,clips,minsize)

K=max(labels);
mapping=zeros(1,K);
k2=0;
for k=1:K
    ct=length(find(labels==k));
    if (ct>=minsize)
        mapping(k)=k2+1; k2=k2+1;
    else
        mapping(k)=0;
    end;
end;
labels2=mapping(labels);
inds=find(labels2>0);
labels2=labels2(inds);
times2=times(inds);
clips2=clips(:,:,inds);

end

function [times2,labels2,clips2]=refine_clusters(times,labels,clips,channels,opts)

[M,T,NC]=size(clips);
M0=length(channels);

for ii=1:3
    K=max(labels);
    WF=zeros(M0,T,K);
    for k=1:K
        inds=find(labels==k);
        WF(:,:,k)=median(clips(channels,:,inds),3);
    end;
    use_it=zeros(1,size(clips,3));
    for k=1:K
        WF0=WF(:,:,k);
        inds=find(labels==k);
        NC0=length(inds);
        clips0=clips(channels,:,inds);
        %clips0_filt=pca_filter_clips(clips0,opts.num_cluster_features);
        
        resid_sumsqr=sum(sum((clips0-repmat(WF0,1,1,NC0)).^2,1),2);
        clips0_sumsqr=sum(sum(clips0.^2,1),2);
        okay_inds=find(resid_sumsqr<clips0_sumsqr);
        if (length(okay_inds)==0) okay_inds=1; end; %just to avoid a possible error. In this case, this cluster will be removed anyway.
        use_it(inds(okay_inds))=1;
    end;
    fprintf(' refine:%d/%d ',sum(use_it),length(use_it));
    inds00=find(use_it);
    times=times(inds00);
    labels=labels(inds00);
    clips=clips(:,:,inds00);
end;

[times,labels,clips]=remove_small_clusters(times,labels,clips,10);

times2=times;
labels2=labels;
clips2=clips;
end

function ret=compute_number_of_templates_to_use(WF,k)
sizes=squeeze(sum(WF.^2,2));
max_sizes=(max(sizes,[],1));
rel_sizes=sizes(k,:)./max_sizes;
labels_to_use=find(rel_sizes>=0.9);
ret=length(labels_to_use);
end

function clips_filt=pca_filter_clips(clips,npca)

[M,T,NC]=size(clips);

[FF,info]=ss_eventfeatures(clips);
FF=FF(1:npca,:);
SS=info.subspace;
SS=SS(:,:,1:npca);

SS0=reshape(SS,[M*T,size(SS,3)]);
clips_filt=reshape(SS0*FF,[M,T,NC]);

end

function labels=do_cluster(clips,opts)
FF=ss_eventfeatures(clips);
FF=FF(1:opts.num_cluster_features,:);
%fprintf('  iso--');
labels0=isosplit(FF);
K=max(labels0);
%fprintf('%d ',K);
labels=labels0;
% labels=zeros(size(labels0));
% lab=0;
% for k=1:K
%     inds=find(labels0==k);
%     fprintf('iso%d--',k);
%     FF2=ss_eventfeatures(clips(:,:,inds));
%     FF2=FF2(1:opts.num_cluster_features,:);
%     labels2=isosplit(FF2);
%     K2=max(labels2);
%     fprintf('%d ',K2);
%     for kk=1:K2
%         inds2=find(labels2==kk);
%         lab=lab+1;
%         labels(inds(inds2))=lab;
%     end;
% end;
end
