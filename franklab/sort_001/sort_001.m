function sort_001(raw_path,output_path,sort_opts)

if nargin<1 test_sort_001; return; end;

def_sort_opts.clip_size=120;
def_sort_opts.shell.min_section_count=200;
def_sort_opts.shell.num_sections_per_shell=4;
def_sort_opts.shell.section_increment=0.5;
def_sort_opts.shell.merge_threshold=0.6;
def_sort_opts.shell.num_features=12;
def_sort_opts.isosplit.isocut_threshold=1.5;
def_sort_opts.filter.samplefreq=30000;
def_sort_opts.filter.freq_min=100;
def_sort_opts.filter.freq_max=10000;
def_sort_opts.filter.outlier_threshold=500;
def_sort_opts.detect.detect_threshold=4;
def_sort_opts.detect.detect_interval=15;
def_sort_opts.plausibility_threshold=2.5;
def_sort_opts.detectibility_threshold=1.8;
def_sort_opts.min_cluster_size=10;

if nargin<3 sort_opts=struct; end;
sort_opts=use_default_opts(sort_opts,def_sort_opts);
sort_opts.shell=use_default_opts(sort_opts.shell,def_sort_opts.shell);
sort_opts.isosplit=use_default_opts(sort_opts.isosplit,def_sort_opts.isosplit);
sort_opts.filter=use_default_opts(sort_opts.filter,def_sort_opts.filter);
sort_opts.detect=use_default_opts(sort_opts.detect,def_sort_opts.detect);

path0=output_path;

%%%% Preprocessing
mscmd_bandpass_filter(raw_path,[path0,'/pre1.mda'],sort_opts.filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],struct);
mscmd_detect([path0,'/pre2.mda'],[path0,'/detect.mda'],sort_opts.detect);
mscmd_extract_clips([path0,'/pre2.mda'],[path0,'/detect.mda'],[path0,'/clips.mda'],struct('clip_size',sort_opts.clip_size));

%%%% Reading clips
fprintf('Reading...\n');
clips=readmda([path0,'/clips.mda']);
[M,T,NC]=size(clips);
clips=clips-repmat(mean(clips,2),1,T,1); %subtract mean over time

%%%% Shell cluster
fprintf('Shell cluster...\n');
sort_opts.shell.isosplit=sort_opts.isosplit;
[labels1,peaks]=shell_cluster(clips,sort_opts.shell);
K=max(labels1);
firings1=zeros(4,NC);
detect=readmda([path0,'/detect.mda']);
times=detect(2,:);
firings1(1:2,:)=detect;
firings1(3,:)=labels1;
firings1(4,:)=peaks;
writemda(firings1,[path0,'/firings1.mda']);

%%%% Split clusters
fprintf('Split clusters...\n');
[firings1_split,cluster_map]=split_clusters_by_peak(clips,firings1);
writemda(firings1_split,[path0,'/firings1_split.mda']);

%%%% Define rclips
fprintf('Define rclips...\n');
o_detect_rclips=sort_opts.detect; o_detect_rclips.detect_threshold=3;
mscmd_detect([path0,'/pre2.mda'],[path0,'/detect_rclips.mda'],o_detect_rclips);
mscmd_extract_clips([path0,'/pre2.mda'],[path0,'/detect_rclips.mda'],[path0,'/rclips.mda'],struct('clip_size',sort_opts.clip_size));
rclips=readmda([path0,'/rclips.mda']);
num_rclips=size(rclips,3);
interval=ceil(num_rclips/5000);
rclips=rclips(:,:,1:interval:end);

%%%% Plausibility scores
fprintf('Plausibility scores...\n');
[firings2_split]=compute_plausibility_scores(clips,firings1_split,rclips);
writemda(firings2_split,[path0,'/firings2_split.mda']);

%%%% Remove events with large plausibility score
firings3_split=firings2_split;
inds0=find(firings2_split(5,:)>sort_opts.plausibility_threshold);
firings3_split(3,inds0)=0;
writemda(firings3_split,[path0,'/firings3_split.mda']);

%%%% Detectibility scores
fprintf('Detectibility scores...\n');
labels3_split=firings3_split(3,:);
templates3_split=ms_templates(clips,labels3_split);
detectibility_scores=compute_detectibility_scores(templates3_split,rclips);

%%%% Remove clusters with low detectibility scores
K_split=length(detectibility_scores);
for k=1:K_split
    if (detectibility_scores(k)<sort_opts.detectibility_threshold)
        firings3_split(3,find(labels3_split==k))=0;
    end;
end;
writemda(firings3_split,[path0,'/firings3_split_B.mda']);

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
writemda(firings,[path0,'/firings.mda']);

fprintf('.\n');

end

function detectibility_scores=compute_detectibility_scores(templates,rclips)
K=size(templates,3);
num_rclips=size(rclips,3);
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
%figure; plot(1:K,template_norms,'b.',1:K,detectibility_scores,'r.','markersize',8);
figure; plot(template_norms,detectibility_scores,'k.','markersize',8);
xlabel('Template norm'); ylabel('Detectibility score');
end

function firings2=compute_plausibility_scores(clips,firings,rclips)
[M,T,L]=size(clips);
[~,~,num_rclips]=size(rclips);
labels=firings(3,:);
K=max(labels);
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

function [labels,clip_peaks]=shell_cluster(clips,opts)
[M,T,NC]=size(clips);

% Compute peaks
fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

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
    else
        CC.labels=[];
        CC.K=0;
    end;
    clusterings{end+1}=CC;
end

% Merge clusters in adjacent shells
for ii=1:length(clusterings)
    K1=clusterings{ii}.K;
    clusterings{ii}.clusters={};
    for k1=1:K1
        clusterings{ii}.clusters{k1}.inds=clusterings{ii}.inds(find(clusterings{ii}.labels==k1));
    end;
    if ii>1
        K2=clusterings{ii-1}.K;
        match_counts=zeros(K1,K2);
        for k1=1:K1
            for k2=1:K2
                inds_k1=clusterings{ii}.inds(find(clusterings{ii}.labels==k1));
                inds_k2=clusterings{ii-1}.inds(find(clusterings{ii-1}.labels==k2));
                match_counts(k1,k2)=length(intersect(inds_k1,inds_k2));
            end;
        end;
        for k1=1:K1
            for k2=1:K2
                numer=match_counts(k1,k2);
                denom=sum(match_counts(k1,:))+sum(match_counts(:,k2))-match_counts(k1,k2);
                if ((denom)&&(numer/denom>=opts.merge_threshold))
                    pct=numer/denom;
                    fprintf('Merging [%d,%d] to [%d,%d] (%d%%)\n',ii,k1,ii-1,k2,floor(pct*100));
                    inds_k1=clusterings{ii}.clusters{k1}.inds;
                    inds_k2=clusterings{ii-1}.clusters{k2}.inds;
                    clusterings{ii}.clusters{k1}.inds=union(inds_k1,inds_k2);
                    clusterings{ii-1}.clusters{k2}.inds=[];
                end;
            end;
        end;
    end;
end;

% Gather clusters
clusters={};
for ii=1:length(clusterings)
    for jj=1:length(clusterings{ii}.clusters)
        CC=clusterings{ii}.clusters{jj};
        if (length(CC.inds)>0)
            clusters{end+1}=CC;
        end;
    end;
end;

% Define labels
labels=zeros(1,NC);
for jj=1:length(clusters)
    labels(clusters{jj}.inds)=jj;
end;

% Redefine labels with no gaps
K=max(labels);
used=zeros(1,K); used(labels)=1; iii=find(used);
mapping=zeros(1,K);
for jj=1:length(iii) mapping(iii(jj))=jj; end;
labels=mapping(labels);
K=max(labels);

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

%sort
[~,inds]=sort(peak_mins);
peak_mins=peak_mins(inds);
peak_maxs=peak_maxs(inds);

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

function opts=use_default_opts(opts,def_opts)
names=fieldnames(def_opts);
for ii=1:length(names)
    if (~isfield(opts,names{ii}))
        opts.(names{ii})=def_opts.(names{ii});
    end;
end;
end

function test_sort_001

close all;

tetrode_num=2;

%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode'];
path0=[mfile_path,sprintf('/output_tetrode%d',tetrode_num)];
if ~exist(path0,'dir') mkdir(path0); end;

%%%% Extract raw data
test_extract_raw_data(raw_path,path0,tetrode_num);

%%%% Sort
sort_001([path0,'/pre0.mda'],path0);

%%%% View output
mv.mode='overview2';
mv.raw=[path0,'/pre2.mda'];
mv.firings=[path0,'/firings.mda'];
ms_mountainview(mv);

end



