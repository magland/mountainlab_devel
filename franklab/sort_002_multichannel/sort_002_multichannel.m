function [firings_path,pre_path]=sort_002_multichannel(raw_path,output_path,sort_opts)
%SORT_001 - Version 001 of sorting based on shell method and isosplit2
%
% Syntax:  firings_path=sort_002_multichannel(raw_path,output_path,sort_opts)
%
% Inputs:
%    raw_path - path to .mda of MxN raw signal data
%    output_path - path to existing DIRECTORY where all output will be written
%    sort_opts - (optional) sorting options, see def_sort_opts in this
%                   script
%
% Outputs:
%    firings_path - path to the firings.mda output file
%    pre_path - path to the preprocessed raw data file
%
% Other m-files required: isosplit2, mscmd_*, ms_*

% Author: Jeremy Magland
% Feb 2016; Last revision: 23-Feb-2016

if nargin<1 test_sort_001_multichannel; return; end;

def_sort_opts.clip_size=120;
def_sort_opts.shell.min_section_count=50;
def_sort_opts.shell.num_sections_per_shell=4;
def_sort_opts.shell.section_increment=0.25;
def_sort_opts.shell.merge_threshold=0.8;
def_sort_opts.shell.num_features=12;
def_sort_opts.isosplit.isocut_threshold=1.5;
def_sort_opts.isosplit.verbose=0;
def_sort_opts.filter.samplefreq=30000;
def_sort_opts.filter.freq_min=100;
def_sort_opts.filter.freq_max=10000;
def_sort_opts.filter.outlier_threshold=500;
def_sort_opts.detect.detect_threshold=4;
def_sort_opts.detect.detect_interval=25;
def_sort_opts.plausibility_threshold=2.5;
def_sort_opts.detectability_threshold=2;
def_sort_opts.min_cluster_size=10;
def_sort_opts.detect.individual_channels=1;
def_sort_opts.adjacency_matrix=[];

if nargin<3 sort_opts=struct; end;
sort_opts=set_default_opts(sort_opts,def_sort_opts);

if (~sort_opts.detect.individual_channels)
    error('individual_channels must be set to 1');
end;

for m=1:size(sort_opts.adjacency_matrix,1)
    sort_opts.adjacency_matrix(m,m)=0;
end;

path0=output_path;

%%%% Preprocessing
mscmd_bandpass_filter(raw_path,[path0,'/pre1.mda'],sort_opts.filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],struct);
mscmd_detect([path0,'/pre2.mda'],[path0,'/detect.mda'],sort_opts.detect);
detect=readmda([path0,'/detect.mda']);

%%%% Sort individual channels
pre2=readmda([path0,'/pre2.mda']);
M=max(detect(1,:));
channels_to_sort=ones(1,M);
if (isfield(sort_opts,'debug_channels'))
    channels_to_sort(:)=0;
    channels_to_sort(sort_opts.debug_channels)=1;
end;
if (isempty(sort_opts.adjacency_matrix))
    sort_opts.adjacency_matrix=ones(M,M);
end;
channel_firings={};
channel_inds={};
for m=1:M
    inds_m=find(detect(1,:)==m);
    detect_m=detect(:,inds_m);
    channel_inds{m}=inds_m;
    fprintf('Sorting channel %d...\n',m);
    neighborhood=[m,find(sort_opts.adjacency_matrix(m,:))];
    if (channels_to_sort(m))
        channel_firings{m}=sort_channel(path0,pre2,detect_m,neighborhood,sort_opts);
    else
        channel_firings{m}=zeros(5,0);
    end;
end;

%%%% Consolidate
fprintf('Consolidating channel firings...\n');
for m=1:M
    channel_firings{m}=consolidate_channel_firings(pre2,channel_firings{m},m,sort_opts);
end;

%%%% Assemble firings
fprintf('Assembling firings...\n');
firings=zeros(5,0);
k_offset=0;
for m=1:M
    firings_m=channel_firings{m};
    labels_m=firings_m(3,:);
    if (~isempty(labels_m))
        K_m=max(labels_m); 
        labels_m=labels_m+k_offset;
        k_offset=k_offset+K_m;
        firings_m(3,:)=labels_m;
        firings=cat(2,firings,firings_m);
    end;
end;
[~,sort_inds]=sort(firings(2,:));
firings=firings(:,sort_inds);

%%%% Writing firings.mda
fprintf('Writing firings.mda...\n');
writemda(firings,[path0,'/firings.mda']);
firings_path=[path0,'/firings.mda'];
pre_path=[path0,'/pre2.mda'];

end

function firings=consolidate_channel_firings(pre,firings,m,sort_opts)
times=firings(2,:);
labels=firings(3,:);
K=max(labels);
clips=ms_extract_clips(pre,times,sort_opts.clip_size);
%templates=compute_geometric_median_templates(clips,labels);
templates=ms_templates(clips,labels);
energies=squeeze(sum(templates.^2,2));
[~,inds]=max(energies,[],1);
labels_to_use=find(inds==m);

events_to_use=zeros(1,length(times));
label_map=zeros(1,K);
kk=1;
for k=labels_to_use
    events_to_use(find(labels==k))=1;
    label_map(k)=kk;
    kk=kk+1;
end;

firings=firings(:,find(events_to_use));
firings(3,:)=label_map(firings(3,:));

end

function firings=sort_channel(path0,pre,detect,neighborhood,sort_opts)
times=detect(2,:);

%%%% Extract clips
fprintf('Extract clips...\n');
clips=ms_extract_clips(pre,times,sort_opts.clip_size);
clips=clips(neighborhood,:,:);
[M,T,NC]=size(clips);
clips=clips-repmat(mean(clips,2),1,T,1); %subtract mean over time

%%%% Branch cluster
fprintf('Branch cluster...\n');
sort_opts.shell.isosplit=sort_opts.isosplit;
%[labels1,peaks]=shell_cluster(clips,sort_opts.shell);
[labels1,peaks]=branch_cluster(clips,sort_opts.shell);
K=max(labels1);
firings1=zeros(4,NC);
firings1(1:2,:)=detect;
firings1(3,:)=labels1;
firings1(4,:)=peaks;
writemda(firings1,[path0,'/firings1.mda']);

% %%%% Shell cluster
% fprintf('Shell cluster...\n');
% sort_opts.shell.isosplit=sort_opts.isosplit;
% [labels1,peaks]=shell_cluster(clips,sort_opts.shell);
% K=max(labels1);
% firings1=zeros(4,NC);
% firings1(1:2,:)=detect;
% firings1(3,:)=labels1;
% firings1(4,:)=peaks;
% writemda(firings1,[path0,'/firings1.mda']);

%%%% Split clusters
fprintf('Split clusters...\n');
[firings1_split,cluster_map]=split_clusters_by_peak(clips,firings1);
writemda(firings1_split,[path0,'/firings1_split.mda']);

%%%% Define rclips
fprintf('Define rclips...\n');
o_detect_rclips.detect_interval=30;
o_detect_rclips.detect_threshold=3;
o_detect_rclips.clip_size=sort_opts.clip_size;
rclips_times=ms_detect(pre(neighborhood,:),o_detect_rclips);
rclips=ms_extract_clips(pre(neighborhood,:),rclips_times,sort_opts.clip_size);
num_rclips=size(rclips,3);
interval=ceil(num_rclips/5000);
rclips=rclips(:,:,1:interval:end);

%%%% Plausibility scores
fprintf('Plausibility scores...\n');
[firings2_split]=compute_plausibility_scores(clips,firings1_split,rclips);
writemda(firings2_split,[path0,'/firings2_split.mda']);

%%%% Remove events with large plausibility score
firings3_split=firings2_split;
% inds0=find(firings2_split(5,:)>sort_opts.plausibility_threshold);
% firings3_split(3,inds0)=0;
% writemda(firings3_split,[path0,'/firings3_split.mda']);

%%%% Detectability scores
fprintf('Detectability scores...\n');
labels3_split=firings3_split(3,:);
templates3_split=ms_templates(clips,labels3_split);
detectability_scores=compute_detectability_scores(templates3_split,rclips);

%%%% Remove clusters with low detectability scores
% K_split=length(detectability_scores);
% for k=1:K_split
%     if (detectability_scores(k)<sort_opts.detectability_threshold)
%         firings3_split(3,find(labels3_split==k))=0;
%     end;
% end;
% writemda(firings3_split,[path0,'/firings3_split_B.mda']);

%%%% Recombine
fprintf('Recombine...\n');
firings3=firings3_split;
inds0=find(firings3(3,:)~=0); firings3(3,inds0)=cluster_map(firings3(3,inds0));
writemda(firings3,[path0,'/firings3.mda']);

%%%% Remove small clusters
fprintf('Remove small clusters...\n');
labels=firings3(3,:);
labels_new=zeros(size(labels));
K=max(labels);
K_new=0;
for k=1:K
    count0=length(find(labels==k));
    if (count0>=sort_opts.min_cluster_size)
        K_new=K_new+1;
        labels_new(find(labels==k))=K_new;
    end;
end;
firings=firings3;
firings(3,:)=labels_new;

fprintf('.\n');

end

function opts=set_default_opts(opts,def_opts)
names=fieldnames(def_opts);
for ii=1:length(names)
    if (isstruct(def_opts.(names{ii})))
        if (~isfield(opts,names{ii}))
            opts.(names{ii})=struct;
        end;
        opts.(names{ii})=set_default_opts(opts.(names{ii}),def_opts.(names{ii}));
    else
        if (~isfield(opts,names{ii}))
            opts.(names{ii})=def_opts.(names{ii});
        end;
    end;
end;
end

function detectability_scores=compute_detectability_scores(templates,rclips)
K=size(templates,3);
num_rclips=size(rclips,3);
detectability_scores=zeros(1,K);
template_norms=zeros(1,K);
for k=1:K
    template_k=templates(:,:,k);
    %figure; ms_view_templates(template0);
    inner_products=squeeze(sum(sum(repmat(template_k,1,1,num_rclips).*rclips,1),2));
    mu=mean(inner_products);
    sigma=sqrt(var(inner_products));
    ip_template0=sum(template_k(:).^2);
    detectability_scores(k)=(ip_template0-mu)/sigma;
    template_norms(k)=sqrt(ip_template0);
end;
%figure; plot(1:K,template_norms,'b.',1:K,detectability_scores,'r.','markersize',8);
figure; plot(template_norms,detectability_scores,'k.','markersize',8);
xlabel('Template norm'); ylabel('Detectability score');
end

function firings2=compute_plausibility_scores(clips,firings,rclips)
[M,T,L]=size(clips);
[~,~,num_rclips]=size(rclips);
labels=firings(3,:);
K=max(labels);
if (isempty(K))||(K==0)
    firings2=firings;
    firings2(5,:)=0;
    return;
end;
templates=compute_geometric_median_templates(clips,labels);

plausibility_scores=zeros(K,L);
for k=1:K
    fprintf('.');
    template=templates(:,:,k);
    if (max(isnan(template(:)))==0)
        weights=get_template_weights(template,5);
        rclips_weighted=rclips.*repmat(weights,1,1,num_rclips);
        [FF_rclips_weighted,subspace]=ms_event_features(rclips_weighted,6);
        for j=1:size(FF_rclips_weighted,1)
            factor=sqrt(var(FF_rclips_weighted(j,:)));
            subspace(:,:,j)=subspace(:,:,j)/factor;
        end
        subspace_vecs=clips_to_vecs(subspace);
        clips_weighted=clips.*repmat(weights,1,1,size(clips,3));
        template_weighted=template.*weights;
        diffs=clips_weighted-repmat(template_weighted,1,1,L);
        FF_rclips_weighted=subspace_vecs'*clips_to_vecs(rclips_weighted);
        FF_diffs=subspace_vecs'*clips_to_vecs(diffs);
        plausibility_scores(k,:)=sqrt(sum(FF_diffs.^2,1));
    else
        plausibility_scores(k,:)=inf;
    end;
end;
fprintf('\n');

[plausibility_scores,labels2]=min(plausibility_scores,[],1);
firings2=firings;
firings2(5,:)=plausibility_scores;

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

function V=clips_to_vecs(clips)
[M,T,K]=size(clips);
V=reshape(clips,M*T,K);
end


function templates=compute_geometric_median_templates(clips,labels)
[M,T,L]=size(clips);
K=max(labels);
templates=zeros(M,T,K);
for k=1:K
    templates(:,:,k)=compute_geometric_median_template(clips(:,:,find(labels==k)));
end;
end

function template=compute_geometric_median_template(clips)
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

function [firings_split,cluster_map]=split_clusters_by_peak(clips,firings)

labels=firings(3,:);
peaks=firings(4,:);

[M,T,L]=size(clips);
K=max(labels);

oo.min_section_count=200;
oo.num_sections_per_shell=1;
oo.section_increment=1;

labels_split=zeros(1,L);
k_split=1;
cluster_map=[];
for k=1:K
    inds_k=find(labels==k);
    [peak_mins,peak_maxs]=define_shells(peaks(inds_k),oo);
    for j=1:length(peak_mins)
        inds0=find((peak_mins(j)<=peaks(inds_k))&(peaks(inds_k)<peak_maxs(j)));
        labels_split(inds_k(inds0))=k_split;
        cluster_map(k_split)=k;
        k_split=k_split+1;
    end;
end;

firings_split=firings;
firings_split(3,:)=labels_split;

end

function [labels,clip_peaks,base_branch_inds]=branch_cluster(clips,opts)
[M,T,NC]=size(clips);

% Compute peaks
fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(1,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(1,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>=abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

%Check for positives and negatives. If so, branch.
inds_neg=find(clip_peaks<0);
inds_pos=find(clip_peaks>=0);
if (length(inds_neg)>0)&&(length(inds_pos)>0)
    [labels_neg]=branch_cluster(clips(:,:,inds_neg),opts);
    [labels_pos]=branch_cluster(clips(:,:,inds_pos),opts);
    labels=zeros(1,NC);
    labels(inds_neg)=labels_neg;
    Kneg=max(1,max(labels_neg));
    labels(inds_pos)=Kneg+labels_pos;
    base_branch_inds=[];
    return;
end;

fprintf('features...\n');
FF=ms_event_features(clips,opts.num_features);
FF=normalize_features(FF);
fprintf('isosplit...\n');
labels0=isosplit2(FF,opts.isosplit);
K0=max(labels0);
fprintf('K0=%d\n',K0);
if (K0>=2)
    %If more than one cluster, then divide and conquer
    labels=zeros(1,NC);
    k_offset=0;
    for k=1:K0
        inds_k=find(labels0==k);
        labels_k=branch_cluster(clips(:,:,inds_k),opts);
        labels(inds_k)=labels_k+k_offset;
        k_offset=k_offset+max(labels_k);
    end;
    base_branch_inds=[];
else %K0=1
    peak_cutoff=0;
    while 1
        if (length(find(abs(clip_peaks)<=peak_cutoff))>=opts.min_section_count)
            if (length(find(abs(clip_peaks)>peak_cutoff))>=opts.min_section_count)
                break;
            end;
        end;
        if (peak_cutoff>max(abs(clip_peaks)))
            break;
        end;
        peak_cutoff=peak_cutoff+opts.section_increment;
    end;
    if (peak_cutoff<max(abs(clip_peaks)))
        fprintf('Using peak cutoff %g\n',peak_cutoff);
        %okay, let's proceed with next section
        inds_low_peak=find(abs(clip_peaks)<=peak_cutoff);
        inds_high_peak=find(abs(clip_peaks)>peak_cutoff);
        [labels_high,~,base_branch_inds_high]=branch_cluster(clips(:,:,inds_high_peak),opts);
        if (length(base_branch_inds_high)>0)
            if (length(find(labels_high(base_branch_inds_high)~=1)))
                error('Unexpected problem. Base branch labels should be 1');
            end;
        end;
        labels=zeros(1,NC);
        base_branch_inds=[inds_low_peak,inds_high_peak(base_branch_inds_high)];
        labels(base_branch_inds)=1;
        labels(inds_high_peak)=labels_high;
    else
        %okay, we are done
        labels=ones(1,NC);
        base_branch_inds=1:NC; %the whole thing is a base
    end;
end;

end

function FF=normalize_features(FF)
norms=sqrt(sum(FF.^2,1));
inds=find(norms~=0);
FF(:,inds)=FF(:,inds)./repmat(norms(inds),size(FF,1),1);
end

function [labels,clip_peaks]=shell_cluster(clips,opts)
[M,T,NC]=size(clips);

% Compute peaks
fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(1,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(1,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>=abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

% Define shells
[peak_mins,peak_maxs]=define_shells(clip_peaks,opts);

% Cluster in each shell
clusterings={};
for ii=1:length(peak_mins)
    rr1=peak_mins(ii);
    rr2=peak_maxs(ii);
    fprintf('rr1=%g rr2=%g... ',rr1,rr2);
    inds_shell=find((clip_peaks>=rr1)&(clip_peaks<rr2));
    CC.inds=inds_shell;
    CC.rr1=rr1;
    CC.rr2=rr2;
    if (length(inds_shell)>1)
        clips_shell=clips(:,:,inds_shell);
        fprintf('features... ');
        FF_shell=ms_event_features(clips_shell,opts.num_features);
        fprintf('isosplit... ');
        labels_shell=isosplit2(FF_shell,opts.isosplit);
        K=max(labels_shell);
        fprintf('K=%d\n',K);
        CC.labels=labels_shell;
        CC.K=K;
        CC.FF_shell=FF_shell;
        CC.clips_shell=clips_shell;
    else
        CC.labels=[];
        CC.K=0;
        CC.FF_shell=zeros(opts.num_features,0);
        CC.clips_shell=zeros(M,T,0);
    end;
    clusterings{end+1}=CC;
end

% Compute all templates and keep track of levels and indices
templates=zeros(M,T,0);
levels=zeros(1,0);
indices=zeros(1,0);
for j=1:length(clusterings)
    CC=clusterings{j};
    fprintf('Computing geometric median templates %d/%d\n',j,length(clusterings));
    templates0=compute_geometric_median_templates(CC.clips_shell,CC.labels);
    K0=size(templates0,3);
    templates=cat(3,templates,templates0);
    levels=cat(2,levels,ones(1,K0)*j);
    indices=cat(2,indices,1:K0);
end;
KK=length(levels);

% Assemble similarity matrix
fprintf('Assemble similarity matrix...\n');
SM=zeros(KK,KK);
for k1=1:KK
    for k2=1:KK
        SM(k1,k2)=compute_template_similarity(templates(:,:,k1),templates(:,:,k2));
    end;
end;
similarity_threshold=0.8;
SM=(SM>=similarity_threshold);
fprintf('downstream transitive fill...\n');
SM=downstream_transitive_fill(SM,levels);
fprintf('Greedy cliques...\n');
CL=ms_greedy_cliques(SM);
labels=zeros(1,NC);
for ii=1:length(CL)
    CL0=CL{ii};
    CL0=sort(CL0); %sort from low to high so that we give preference to higher peak classification
    for jj=1:length(CL0)
        aa=CL0(jj);
        level0=levels(aa);
        index0=indices(aa);
        CC=clusterings{level0};
        inds0=find(CC.labels==index0);
        labels(CC.inds(inds0))=ii;
    end;
end;

figure; ms_view_templates_from_clips(clips,labels);
drawnow;

end

function ret=compute_template_similarity(T1,T2)
T1=T1-repmat(mean(T1,2),1,size(T1,2));
T2=T2-repmat(mean(T2,2),1,size(T2,2));
T1=reshape(T1,1,length(T1(:)));
T2=reshape(T2,1,length(T2(:)));
norm1=sqrt(T1*T1');
norm2=sqrt(T2*T2');
T1=T1/norm1;
T2=T2/norm2;
ret=T1*T2';
end

function A=downstream_transitive_fill(A,levels)
K=length(levels);
something_changed=1;
while something_changed
something_changed=0;
for k1=1:K
    for k2=1:K
        for k3=1:K
            if (A(k1,k2)&&(A(k2,k3))&&(~A(k1,k3)))
                if ((levels(k1)<levels(k2))&&(levels(k2)<levels(k3)))
                    A(k1,k3)=1;
                    A(k3,k1)=1;
                    something_changed=1;
                end;
            end;
        end;
    end;
end;
end;
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
        count0=length(find((min0<=clip_peaks)&(clip_peaks<max0)));
        if (count0>=opts.min_section_count)
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
                count0=length(find((min0<=clip_peaks)&(clip_peaks<max0)));
                if (count0>0)
                    section_mins_neg=[section_mins_neg,-inf];
                    section_maxs_neg=[section_maxs_neg,max0];
                end;
                max0=-inf; min0=-inf;
            end;
        end;
    end;
end;

% do the positives
section_mins_pos=[];
section_maxs_pos=[];
if (max(clip_peaks)>0)
    min0=0;
    max0=min0+opts.section_increment;
    while 1
        if (min0>max(clip_peaks)) break; end;
        count0=length(find((min0<=clip_peaks)&(clip_peaks<max0)));
        if (count0>=opts.min_section_count)
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
                count0=length(find((min0<=clip_peaks)&(clip_peaks<max0)));
                if (count0>0)
                    section_mins_pos=[section_mins_pos,min0];
                    section_maxs_pos=[section_maxs_pos,inf];
                end;
                min0=inf; max0=inf;
            end;
        end;
    end;
end;

peak_mins=[];
peak_maxs=[];
%the negatives
for j=1:length(section_mins_neg)
    peak_maxs(end+1)=section_maxs_neg(j);
    ii=j+opts.num_sections_per_shell-1;
    if (ii<=length(section_mins_neg))
        peak_mins(end+1)=section_mins_neg(ii);
    else
        peak_mins(end+1)=-inf;
    end;
end;
%the positives
for j=1:length(section_mins_pos)
    peak_mins(end+1)=section_mins_pos(j);
    ii=j+opts.num_sections_per_shell-1;
    if (ii<=length(section_maxs_pos))
        peak_maxs(end+1)=section_maxs_pos(ii);
    else
        peak_maxs(end+1)=inf;
    end;
end;


%sort (actually, don't sort!)
%[~,inds]=sort(peak_mins);
%peak_mins=peak_mins(inds);
%peak_maxs=peak_maxs(inds);

end

function test_extract_raw_data(raw_path,output_path,tetrode_num)

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


function test_sort_001_multichannel

close all;

tetrode_num=1;

%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode'];
path0=[mfile_path,sprintf('/output_tetrode%d',tetrode_num)];
if ~exist(path0,'dir') mkdir(path0); end;

%%%% Extract raw data
test_extract_raw_data(raw_path,path0,tetrode_num);

%%%% Sort
sort_001_multichannel([path0,'/pre0.mda'],path0);

%%%% View output
mv.mode='overview2';
mv.raw=[path0,'/pre2.mda'];
mv.firings=[path0,'/firings.mda'];
ms_mountainview(mv);

end


