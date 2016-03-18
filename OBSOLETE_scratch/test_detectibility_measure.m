function test_detectibility_measure

close all;

clip_size=100;

% Load data
path0='franklab/test_hippocampal_02_10_2016/output_tetrode1';
pre2=readmda([path0,'/pre2.mda']);
firings=readmda([path0,'/firings.mda']);
times=firings(2,:);
labels=firings(3,:);
peaks=firings(4,:);
%inds=find((abs(peaks)<6)&(abs(peaks)>5));
inds=find(abs(peaks)>0);
times=times(inds); labels=labels(inds); peaks=peaks(inds);
labels=split_clusters_by_peak_amplitudes(peaks,labels,struct('section_increment',1,'min_section_count',200,'num_sections_per_shell',1));
clips=ms_extract_clips(pre2,times,clip_size);
[M,T,L]=size(clips);
clips=clips-repmat(mean(clips,2),1,T,1);
%templates=ms_templates(clips,labels);
templates=geometric_median_templates(clips,labels);
figure; ms_view_templates(templates);
K=size(templates,3);

%Define rclips
ttt=ms_detect(pre2,struct('detect_threshold',3,'detect_interval',15,'clip_size',clip_size));
interval=ceil(length(ttt)/5000);
ttt=ttt(1:interval:end);
rclips=ms_extract_clips(pre2,ttt,clip_size);
num_rclips=size(rclips,3);
fprintf('%d reference clips\n',num_rclips);

% Compute detectibility scores
detectibility_scores=zeros(1,K);
template_norms=zeros(1,K);
for k=1:K
    template_k=templates(:,:,k);
    %figure; ms_view_templates(template0);
    inner_products=squeeze(sum(sum(repmat(template_k,1,1,num_rclips).*rclips,1),2));
    mu=mean(inner_products);
    sigma=sqrt(var(inner_products));
    ip_template0=sum(template_k(:).^2);
    detectibility_scores(k)=(ip_template0-mu)/sigma;
    template_norms(k)=sqrt(ip_template0);
end;

% Compute plausibility scores
plausibility_scores=zeros(K,L);
for k=1:K
    if (~isnan(detectibility_scores(k)))
        fprintf('.');
        template=templates(:,:,k);
        weights=get_template_weights(template,5);
        rclips_weighted=rclips.*repmat(weights,1,1,num_rclips);
        [FF_rclips_weighted,subspace]=ms_event_features(rclips_weighted,6);
        for j=1:size(FF_rclips_weighted,1)
            factor=sqrt(var(FF_rclips_weighted(j,:)));
            subspace(:,:,j)=subspace(:,:,j)/factor;
        end
        subspace_vecs=templates_to_vecs(subspace);
        clips_weighted=clips.*repmat(weights,1,1,size(clips,3));
        template_weighted=template.*weights;
        diffs=clips_weighted-repmat(template_weighted,1,1,L);
        FF_rclips_weighted=subspace_vecs'*templates_to_vecs(rclips_weighted);
        FF_diffs=subspace_vecs'*templates_to_vecs(diffs);
        plausibility_scores(k,:)=sqrt(sum(FF_diffs.^2,1));
    else
        plausibility_scores(k,:)=inf;
    end;
    %figure; ms_view_clusters(cat(2,FF_rclips_weighted,FF_diffs),cat(2,ones(1,num_rclips),ones(1,L)*2));
%     inds=find(labels==k);
%     inds_not=find(labels~=k);
%     %figure; ms_view_clusters(cat(2,FF_diffs(:,inds_not),FF_diffs(:,inds)),cat(2,ones(1,length(inds_not)),ones(1,length(inds))*2));
%     figure; ms_view_clusters(FF_diffs(:,inds));
%     hold on; [xx,yy,zz]=sphere; h=surfl(xx,yy,zz); set(h, 'FaceAlpha', 0.1);
%     title(sprintf('k=%d, detectibility=%g',k,detectibility_scores(k)));
end;
fprintf('\n');

detectibility_threshold=2;
k0=1;
for k=1:K
    if (detectibility_scores(k)>=detectibility_threshold)&&(~isnan(detectibility_scores(k)))
        labels(find(labels==k))=k0;
        k0=k0+1;
    else
        labels(find(labels==k))=0;
    end;
end;

templates=geometric_median_templates(clips,labels);
figure; ms_view_templates(templates);

figure; plot(1:K,template_norms,'b.',1:K,detectibility_scores,'r.','markersize',8);
figure; plot(template_norms,detectibility_scores,'k.','markersize',8);

mv.raw=pre2;
mv.firings=TL2F(times,labels);
ms_mountainview(mv);

end

function V=templates_to_vecs(templates)
[M,T,K]=size(templates);
V=reshape(templates,M*T,K);
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

function [peak_mins,peak_maxs]=define_shells(clip_peaks,opts)

% do the negatives
section_mins_neg=[];
section_maxs_neg=[];
if (min(clip_peaks)<0)
    max0=0;
    min0=max0-opts.section_increment;
    while 1
        if (max0<min(clip_peaks)) break; end;
        if (length(find((min0<=clip_peaks)&(clip_peaks<max0)))>=opts.min_section_count)
            if (length(find(clip_peaks<min0))>=opts.min_section_count)
                section_mins_neg=[section_mins_neg,min0];
                section_maxs_neg=[section_maxs_neg,max0];
                max0=min0;
                min0=max0-opts.section_increment;
            else
                section_mins_neg=[section_mins_neg,-inf];
                section_maxs_neg=[section_maxs_neg,max0];
                max0=-inf; min0=-inf;
            end;
        else
            min0=min0-opts.section_increment;
            if (min0<min(clip_peaks))
                section_mins_neg=[section_mins_neg,-inf];
                section_maxs_neg=[section_maxs_neg,max0];
                max0=-inf; min0=-inf;
            end;
        end;
    end;
end;
% do the overlap
for j=1:length(section_maxs_neg)-1
    section_mins_neg(j)=section_mins_neg(j+1);
end;

% do the positives
section_mins_pos=[];
section_maxs_pos=[];
if (max(clip_peaks)>0)
    min0=0;
    max0=min0+opts.section_increment;
    while 1
        if (min0>max(clip_peaks)) break; end;
        if (length(find((min0<=clip_peaks)&(clip_peaks<max0)))>=opts.min_section_count)
            if (length(find(clip_peaks>=max0))>=opts.min_section_count)
                section_mins_pos=[section_mins_pos,min0];
                section_maxs_pos=[section_maxs_pos,max0];
                min0=max0;
                max0=min0+opts.section_increment;
            else
                section_mins_pos=[section_mins_pos,min0];
                section_maxs_pos=[section_maxs_pos,inf];
                min0=inf; max0=inf;
            end;
        else
            max0=max0+opts.section_increment;
            if (max0>max(clip_peaks))
                section_mins_pos=[section_mins_pos,min0];
                section_maxs_pos=[section_maxs_pos,inf];
                min0=inf; max0=inf;
            end;
        end;
    end;
end;
% do the overlap
for j=1:length(section_mins_pos)-1
    section_maxs_pos(j)=section_maxs_pos(j+1);
end;

% combine and sort
section_mins=[section_mins_neg,section_mins_pos];
section_maxs=[section_maxs_neg,section_maxs_pos];
[~,inds]=sort(section_mins);
section_mins=section_mins(inds);
section_maxs=section_maxs(inds);

peak_mins=section_mins;
for j=2:opts.num_sections_per_shell
    peak_mins=[-inf,peak_mins];
end;
peak_maxs=section_maxs;
for j=2:opts.num_sections_per_shell
    peak_maxs=[peak_maxs,inf];
end;

end

function [labels_out]=split_clusters_by_peak_amplitudes(clip_peaks,labels,opts)
NC=length(clip_peaks);
K=max(labels);
labels_out=zeros(1,NC);
for k=1:K
    inds=find(labels==k);
    labels0=split_cluster_by_peak_amplitudes(clip_peaks(inds),opts);
    labels_out(inds)=labels0+max(labels_out);
end;
end

function [labels]=split_cluster_by_peak_amplitudes(clip_peaks,opts)
NC=length(clip_peaks);
[peak_mins,peak_maxs]=define_shells(clip_peaks,opts);

labels=zeros(1,NC);
for ii=1:length(peak_mins)
    inds=find((clip_peaks>=peak_mins(ii))&(clip_peaks<peak_maxs(ii)));
    labels(inds)=max(labels)+1;
end;

end

function templates=geometric_median_templates(clips,labels)
K=max(labels);
[M,T,NC]=size(clips);
templates=zeros(M,T,K);
for k=1:K
    inds_kk=find(labels==k);
    clips_kk=clips(:,:,inds_kk);
    templates(:,:,k)=compute_clips_template(clips_kk);
end;
end

function template=compute_clips_template(clips)
[M,T,NC]=size(clips);
if (length(clips(:))==0)
    template=zeros(M,T);
    return;
end;
[M,T,NC]=size(clips);
num_features=18;
FF=ms_event_features(clips,num_features);
FFmm=ms_geometric_median(FF);
diffs=FF-repmat(FFmm,1,NC);
dists=sqrt(sum(diffs.^2,1));
sorted_dists=sort(dists);
dist_cutoff=sorted_dists(ceil(length(sorted_dists)*0.3));
inds=find(dists<=dist_cutoff);
template=mean(clips(:,:,inds),3);
end

