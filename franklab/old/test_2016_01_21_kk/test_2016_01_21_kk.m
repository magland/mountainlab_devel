function test_2016_01_21_kk

close all;

mfile_path=fileparts(mfilename('fullpath'));
output_kk_path=sprintf('%s/output_kk',mfile_path);
output_ms_path=sprintf('%s/../experiment1/output',mfile_path);
if (~exist(output_kk_path,'dir')) mkdir(output_kk_path); end;
kwik_path=sprintf('%s/../raw/2016_01_21/ms11d45_new.kwik',mfile_path);
locations_path=sprintf('%s/locations.mda',output_ms_path);

channel_subset=[136,141,147,142,145,107,148,150,152,153,98,99,24,154,19,100,8,25,70,127,155,156,133,138,79,17,14,12,157];

opts=get_sorting_options;

% Read the KlustaKwik data, if not already done
if ~exist(sprintf('%s/clusters.mda',output_kk_path),'file')
    fprintf('Read KlustaKwik...\n');
    times=double(hdf5read(kwik_path,'/channel_groups/2/spikes/time_samples'));
    labels=double(hdf5read(kwik_path,'/channel_groups/2/spikes/clusters/main'));
    inds=find(times<opts.o_extract.t2); %truncate
    times=times(inds);
    labels=labels(inds);
    [~,sort_inds]=sort(times);
    times=times(sort_inds); labels=labels(sort_inds);
    CC_kk=zeros(3,length(times));
    CC_kk(2,:)=times;
    CC_kk(3,:)=labels;
    writemda(CC_kk,sprintf('%s/clusters.mda',output_kk_path));
    
    NE=length(times);
    new_labels=zeros(1,NE);
    for j=1:length(channel_subset)
        inds=find(labels==channel_subset(j));
        new_labels(inds)=j;
    end;
    CC_kk(3,:)=new_labels;
    inds_to_use=find(new_labels>0);
    CC_kk=CC_kk(:,inds_to_use);
    writemda(CC_kk,sprintf('%s/clusters_subset.mda',output_kk_path));
end;

mscmd_templates([output_ms_path,'/pre1.mda'],[output_kk_path,'/clusters.mda'],[output_kk_path,'/templates.mda'],opts.o_templates);
mscmd_extract_clips([output_ms_path,'/pre1.mda'],[output_kk_path,'/clusters.mda'],[output_kk_path,'/clips.mda'],[output_kk_path,'/clips_index.mda'],opts.o_extract_clips);

mscmd_templates([output_ms_path,'/pre1.mda'],[output_kk_path,'/clusters_subset.mda'],[output_kk_path,'/templates_subset.mda'],opts.o_templates);
mscmd_extract_clips([output_ms_path,'/pre1.mda'],[output_kk_path,'/clusters_subset.mda'],[output_kk_path,'/clips_subset.mda'],[output_kk_path,'/clips_subset_index.mda'],opts.o_extract_clips);
mscmd_cross_correlograms([output_kk_path,'/clusters_subset.mda'],[output_kk_path,'/cross_correlograms_subset.mda'],opts.o_cross_correlograms.max_dt);

mscmd_templates([output_ms_path,'/pre1.mda'],[output_ms_path,'/clusters.mda'],[output_kk_path,'/templates_ms.mda'],opts.o_templates);
mscmd_extract_clips([output_ms_path,'/pre1.mda'],[output_ms_path,'/clusters.mda'],[output_kk_path,'/clips_ms.mda'],[output_kk_path,'/clips_ms_index.mda'],opts.o_extract_clips);
mscmd_cross_correlograms([output_ms_path,'/clusters.mda'],[output_kk_path,'/cross_correlograms_ms.mda'],opts.o_cross_correlograms.max_dt);

mscmd_confusion_matrix([output_kk_path,'/clusters_subset.mda'],[output_ms_path,'/clusters.mda'],[output_kk_path,'/confusion_matrix.mda'],3);

CM=readmda([output_kk_path,'/confusion_matrix.mda']);

[K1,K2]=size(CM); K1=K1-1; K2=K2-1;
CM_row_normalized=CM./repmat(sum(CM,2),1,K2+1);
CM_column_normalized=CM./repmat(sum(CM,1),K1+1,1);
CM_normalized=(CM_row_normalized+CM_column_normalized)/2;

reordering=get_best_reordering(CM_normalized');

if 0
    clusters_ms=readmda([output_ms_path,'/clusters.mda']);
    clusters_ms_re=reorder_clusters(clusters_ms,reordering);
    writemda(clusters_ms_re,[output_kk_path,'/clusters_ms_re.mda']);
end;

mscmd_templates([output_ms_path,'/pre1.mda'],[output_kk_path,'/clusters_ms_re.mda'],[output_kk_path,'/templates_ms_re.mda'],opts.o_templates);
mscmd_extract_clips([output_ms_path,'/pre1.mda'],[output_kk_path,'/clusters_ms_re.mda'],[output_kk_path,'/clips_ms_re.mda'],[output_kk_path,'/clips_ms_re_index.mda'],opts.o_extract_clips);
mscmd_cross_correlograms([output_kk_path,'/clusters_ms_re.mda'],[output_kk_path,'/cross_correlograms_ms_re.mda'],opts.o_cross_correlograms.max_dt);


if 0
    clusters=readmda([output_kk_path,'/clusters_subset.mda']);
    templates=readmda([output_kk_path,'/templates_subset.mda']);
    clips=readmda([output_kk_path,'/clips_subset.mda']);
    clips_index=readmda([output_kk_path,'/clips_subset_index.mda']);
    labels=clusters(3,:);

    figure; ms_view_templates(templates);
    %figure; ms_view_templates_from_clips(clips,clips_index,struct('show_stdev',1));
end;

if 0
    clusters=readmda([output_kk_path,'/clusters_ms_re.mda']);
    templates=readmda([output_kk_path,'/templates_ms_re.mda']);
    clips=readmda([output_kk_path,'/clips_ms_re.mda']);
    clips_index=readmda([output_kk_path,'/clips_ms_re_index.mda']);
    labels=clusters(3,:);

    figure; ms_view_templates(templates);
    %figure; ms_view_templates_from_clips(clips,clips_index,struct('show_stdev',1));
end;

if 0
    open_mountainview([output_ms_path,'/pre1.mda'],[output_kk_path,'/clusters_ms_re.mda'],[output_kk_path,'/templates_ms_re.mda'],[output_kk_path,'/clips_ms_re.mda'],[output_kk_path,'/clips_ms_re_index.mda'],[output_kk_path,'/cross_correlograms_ms_re.mda'],locations_path);
end;

if 0
    open_mountainview([output_ms_path,'/pre1.mda'],[output_kk_path,'/clusters_subset.mda'],[output_kk_path,'/templates_subset.mda'],[output_kk_path,'/clips_subset.mda'],[output_kk_path,'/clips_subset_index.mda'],[output_kk_path,'/cross_correlograms_subset.mda'],locations_path);
end;

if 1
    figure; imagesc(CM');
    figure; imagesc(CM_row_normalized(:,[reordering,K2+1])');
    figure; imagesc(CM_column_normalized(:,[reordering,K2+1])');
    figure; imagesc(CM_normalized(:,[reordering,K2+1])');
end;

if 0
    open_mv_compare([output_ms_path,'/pre3.mda'],[output_kk_path,'/clusters_subset.mda'],[output_kk_path,'/clusters_ms_re.mda']);
end

end

function clusters2=reorder_clusters(clusters,reorder)
clusters2=clusters;
[~,reorder2]=sort(reorder);
clusters2(3,:)=reorder2(clusters(3,:));
end

function open_mv_compare(raw,clusters1,clusters2)
mfile_path=fileparts(mfilename('fullpath'));
exe_fname=sprintf('%s/../../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s --mode=compare_labels ',exe_fname)];

cmd=[cmd,sprintf('--raw=%s ',raw)];
cmd=[cmd,sprintf('--cluster=%s ',clusters1)];
cmd=[cmd,sprintf('--cluster2=%s ',clusters2)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));
end

function open_mountainview(raw,clusters,templates,clips,clips_index,cross_correlograms,locations)

mfile_path=fileparts(mfilename('fullpath'));
path0=sprintf('%s/output',mfile_path);
exe_fname=sprintf('%s/../../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s ',exe_fname)];

cmd=[cmd,sprintf('--raw=%s ',raw)];
cmd=[cmd,sprintf('--templates=%s ',templates)];
cmd=[cmd,sprintf('--clips=%s ',clips)];

cmd=[cmd,sprintf('--cluster=%s ',clusters)];
cmd=[cmd,sprintf('--clips-index=%s ',clips_index)];

cmd=[cmd,sprintf('--cross-correlograms=%s ',cross_correlograms)];

cmd=[cmd,sprintf('--locations=%s ',locations)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));

end

function reordering=get_best_reordering_old(templates1,templates2)
[M1,T1,K1]=size(templates1);
[M2,T2,K2]=size(templates2);
reordering=zeros(1,K1);

for j=1:min(K1,K2)
    template_2=templates2(:,:,j);
    ind0=find_best_template_match(templates1,template_2);
    reordering(ind0)=j;
    templates1(:,:,ind0)=0;
end;

unused_inds=find(reordering==0);
reordering(unused_inds)=min(K1,K2)+1:K1;

[~,reordering]=sort(reordering);

end

function reordering=get_best_reordering(CM)
[K1,K2]=size(CM);
K1=K1-1;
K2=K2-1;
HH=Hungarian(-CM(1:K1,1:K2));

map12=zeros(1,K1);
for j=1:K1
    ind0=find(HH(j,1:K2)==1);
    if (length(ind0)==0) map12(j)=0;
    else map12(j)=ind0(1);
    end;
end;

reordering=map12;

used=zeros(1,K1);
used(reordering(reordering>0))=1;
unused_inds=find(used==0);
inds=find(reordering==0);
reordering(inds)=unused_inds;

[~,reordering]=sort(reordering);


end

function ind=find_best_template_match(templates1,template_2)
K=size(templates1,3);
sumsqrs=squeeze(sum(sum(templates1.^2,1),2));
diffs=templates1-repmat(template_2,1,1,K);
diffs=squeeze(sum(sum(diffs.^2,1),2));
inds0=find(sumsqrs==0);
diffs(inds0)=inf;
[~,ind]=min(diffs);
ind=ind(1);
end

function opts=get_sorting_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Raw extract options
opts.o_extract.num_channels=72;
opts.o_extract.channels=[37:52,68,69];
opts.o_extract.t1=0; opts.o_extract.t2=19e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial bandpass filter options
opts.o_filter.samplefreq=30000;
opts.o_filter.freq_min=600;
opts.o_filter.freq_max=4000;
opts.o_filter.outlier_threshold=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Whitening options
opts.o_whiten.ncomp=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-whitening bandpass filter options
opts.o_filter2.samplefreq=30000;
opts.o_filter2.freq_min=600;
opts.o_filter2.freq_max=4000;
opts.o_filter2.outlier_threshold=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detection options
opts.o_detect.inner_window_width=15;
opts.o_detect.outer_window_width=100000;
opts.o_detect.threshold=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature extraction options
opts.o_features.num_features=6;
opts.o_features.clip_size=60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clustering options
opts.o_cluster=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster splitting options
opts.o_split_clusters.num_features=3; 
opts.o_split_clusters.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Template extraction options
opts.o_templates.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster consolidation options
opts.o_consolidate.coincidence_threshold=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit options
opts.o_fit=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross-correlogram options
opts.o_cross_correlograms.max_dt=1500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clip extraction options
opts.o_extract_clips.clip_size=opts.o_templates.clip_size;
end