function test_hippocampal_02_10_2016

%%%% Parameters and settings
tetrode_num=1;
plausibility_threshold=0.6;
merge_threshold=0.8; %Keep this high for now
tt_range=[4,15];
num_tt_steps=12;
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
[labels,peaks]=shell_cluster(clips,tt_range,num_tt_steps,tt_overlap,num_features,merge_threshold);
K=max(labels);
clusters=zeros(4,NC);
detect=readmda([path0,'/detect.mda']);
clusters(1:2,:)=detect;
clusters(3,:)=labels;
clusters(4,:)=peaks;
writemda(clusters,[path0,'/clusters.mda']);
mscmd_templates([path0,'/pre2.mda'],[path0,'/clusters.mda'],[path0,'/templates.mda'],o_templates);
templates=readmda([path0,'/templates.mda']);
figure; 
ms_view_templates(templates);

end

function template=compute_clips_medoid(clips)
num_features=18;
FF=ms_event_features(clips,num_features);
[M,N]=size(FF);
dists=zeros(N,N);
for m=1:M
    [grid1,grid2]=ndgrid(FF(m,:),FF(m,:));
    dists=dists+sqrt((grid1-grid2).^2);
end;
avg_dists=mean(dists,1);
[~,ind]=min(avg_dists);
m=ind(1);

sorted_dists=sort(dists(m,:));
dist_cutoff=sorted_dists(ceil(length(sorted_dists)*0.3));
inds=find(dists(m,:)<=dist_cutoff);
template=mean(clips(:,:,inds),3);

end

function [labels,clip_peaks]=shell_cluster(clips,tt_range,num_tt_steps,tt_overlap,num_features,merge_threshold)
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
    else
        CC.labels=[];
        CC.K=0;
    end;
    clusterings{end+1}=CC;
end

for ii=1:length(clusterings)
    K1=clusterings{ii}.K;
    clusterings{ii}.clusters={};
    labels1=clusterings{ii}.labels;
    for k1=1:K1
        inds_k1=find(labels1==k1);
        CL.inds=clusterings{ii}.inds(inds_k1);
        if (ii>1)
            K2=clusterings{ii-1}.K;
            labels2=clusterings{ii-1}.labels;
            [inds_intersect,ii1,ii2]=intersect(clusterings{ii}.inds,clusterings{ii-1}.inds);
            label_counts=zeros(1,K2);
            label_tots=zeros(1,K2);
            for k2=1:K2
                label_counts(k2)=length(find((labels1(ii1)==k1)&(labels2(ii2)==k2)));
                label_tots(k2)=length(find((labels1(ii1)==k1)|(labels2(ii2)==k2)));
            end;
            inds00=find(label_counts>=label_tots*merge_threshold);
            if (length(inds00)>0)
                k2=inds00(1);
                fprintf('Merging [%d,%d] to [%d,%d]\n',ii,k1,ii-1,k2);
                CL.inds=union(CL.inds,clusterings{ii-1}.inds(find(labels2==k2)));
                clusterings{ii-1}.clusters{k2}.inds=[];
            end;
        end;
        clusterings{ii}.clusters{end+1}=CL;
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

