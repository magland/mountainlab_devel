function test_hippocampal_02_10_2016

close all; drawnow;

%%%% Parameters and settings
tetrode_num=1;
plausibility_threshold=0.6;
merge_threshold=0.8;
num_tt_steps=12;
tt_overlap=1;
num_features=12;
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
o_detect.threshold=4;
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=15;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=60;
o_whiten=struct;
o_templates.clip_size=o_extract_clips.clip_size;

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
[labels,peaks]=shell_cluster(clips,num_tt_steps,tt_overlap,num_features,merge_threshold);
K=max(labels);
clusters=zeros(4,NC);
detect=readmda([path0,'/detect.mda']);
times=detect(2,:);
clusters(1:2,:)=detect;
clusters(3,:)=labels;
clusters(4,:)=peaks;
writemda(clusters,[path0,'/clusters.mda']);
%mscmd_templates([path0,'/pre2.mda'],[path0,'/clusters.mda'],[path0,'/templates.mda'],o_templates);
%templates=readmda([path0,'/templates.mda']);
%figure; 
%ms_view_templates(templates);

%%%% Writing output and preparing view
fprintf('Writing output and preparing view...\n');
pre2=readmda([path0,'/pre2.mda']);
[clips1,clips1_index]=ms_create_clips_index(ms_extract_clips(pre2,times,o_extract_clips.clip_size),labels);
writemda(clips1,[path0,'/clips0.mda']);
writemda(clips1_index,[path0,'/clips0_index.mda']);
writemda(clusters,[path0,'/clusters.mda']);
%writemda(corr_matrix,[path0,'/correlation_matrix.mda']);

%%%% Cross correlograms and templates
mscmd_cross_correlograms([path0,'/clusters.mda'],[path0,'/cross_correlograms.mda'],cross_correlograms_max_dt);
mscmd_templates([path0,'/pre0_mild.mda'],[path0,'/clusters.mda'],[path0,'/templates_raw.mda'],struct('clip_size',200));
mscmd_templates([path0,'/pre2.mda'],[path0,'/clusters.mda'],[path0,'/templates.mda'],struct('clip_size',200));
templates=readmda([path0,'/templates.mda']);
figure; ms_view_templates(templates);

%%%% MountainView
view_params.raw=[path0,'/pre2.mda'];
view_params.clusters=[path0,'/clusters.mda'];
view_params.cross_correlograms=[path0,'/cross_correlograms.mda'];
view_params.templates=[path0,'/templates.mda'];
view_params.clips=[path0,'/clips0.mda'];
view_params.clips_index=[path0,'/clips0_index.mda'];
ms_mountainview(view_params);

%%%% Split clusters by peak amplitudes
templates_split=split_clusters_by_peak_amplitudes(clips,labels);
figure; ms_view_templates(templates_split);

end

function template=compute_clips_template(clips)
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

function [labels,clip_peaks]=shell_cluster(clips,num_tt_steps,tt_overlap,num_features,merge_threshold)
[M,T,NC]=size(clips);

fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

sorted_abs_peaks=sort(abs(clip_peaks));
tt1_list=zeros(1,num_tt_steps);
tt2_list=zeros(1,num_tt_steps);
for jj=1:num_tt_steps
    ind1=max(1,ceil(length(sorted_abs_peaks)*(jj-1)/num_tt_steps+1));
    ind2=min(length(sorted_abs_peaks),ceil(length(sorted_abs_peaks)*(jj-1+1+tt_overlap)/num_tt_steps+1));
    tt1_list(jj)=sorted_abs_peaks(ind1);
    tt2_list(jj)=sorted_abs_peaks(ind2);
end;

clusterings={};
for ii=1:length(tt1_list)
    tt1=tt1_list(ii);
    tt2=tt2_list(ii);
    if (ii==length(tt1_list)), tt2=inf; end;
    fprintf('tt1=%g tt2=%g... ',tt1,tt2);
    inds_tt=find((abs(clip_peaks)>=tt1)&(abs(clip_peaks)<tt2));
    CC.inds=inds_tt;
    if (length(inds_tt)>1)
        clips_tt=clips(:,:,inds_tt);
        fprintf('features... ');
        [FF_tt,subspace_tt]=ms_event_features(clips_tt,num_features);
        fprintf('isosplit... ');
        labels_tt=isosplit2(FF_tt,struct('verbose',0,'verbose3',0));
        K=max(labels_tt);
        fprintf('K=%d\n',K);
        CC.labels=labels_tt;
        CC.K=K;
    else
        CC.labels=[];
        CC.K=0;
    end;
    clusterings{end+1}=CC;
end

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
                if ((denom)&&(numer/denom>=merge_threshold))
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

clusters={};
for ii=1:length(clusterings)
    for jj=1:length(clusterings{ii}.clusters)
        CC=clusterings{ii}.clusters{jj};
        if (length(CC.inds)>0)
            clusters{end+1}=CC;
        end;
    end;
end;

labels=zeros(1,NC);
for jj=1:length(clusters)
    labels(clusters{jj}.inds)=jj;
end;

K=max(labels);
used=zeros(1,K); used(labels)=1; iii=find(used);
mapping=zeros(1,K);
for jj=1:length(iii) mapping(iii(jj))=jj; end;
labels=mapping(labels);
K=max(labels);

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

function templates=split_clusters_by_peak_amplitudes(clips,labels)
[M,T,NC]=size(clips);
templates=zeros(M,T,0);
K=max(labels);
for k=1:K
    inds=find(labels==k);
    templates0=split_cluster_by_peak_amplitudes(clips(:,:,inds));
    templates=cat(3,templates,templates0);
    if (k<K) templates=cat(3,templates,zeros(M,T,2)); end;
end;
end

function [templates,labels]=split_cluster_by_peak_amplitudes(clips)
[M,T,NC]=size(clips);
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

incr=1;

candidate_cutoffs=ceil(min(clip_peaks)/incr)*incr:incr:floor(max(clip_peaks)/incr)*incr;
cutoffs=[-inf];
for ii=1:length(candidate_cutoffs)
    cc=candidate_cutoffs(ii);
    n1=length(find((clip_peaks<cc)&(clip_peaks>=cutoffs(end))));
    n2=length(find(clip_peaks>=cc));
    if (n1>=30)&&(n2>=30)
        cutoffs=[cutoffs,cc];
    end;
end;
cutoffs=[cutoffs,inf];

templates=zeros(M,T,length(cutoffs)-1);
labels=zeros(1,NC);
for ii=1:length(cutoffs)-1
    inds=find((clip_peaks>=cutoffs(ii))&(clip_peaks<cutoffs(ii+1)));
    templates(:,:,ii)=compute_clips_template(clips(:,:,inds));
    labels(inds)=max(labels)+1;
end;

end
