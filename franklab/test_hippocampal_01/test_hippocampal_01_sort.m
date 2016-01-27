function test_hippocampal_01_sort

mfile_path=fileparts(mfilename('fullpath'));
path_raw=sprintf('%s/../raw/hippocampal/tetrode',mfile_path);
path0=sprintf('%s/output',mfile_path);
if ~exist(path0,'dir')
    mkdir(path0);
end;

opts=get_sorting_options;
opts_pre=get_preprocessing_options;

raw_mat_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mat',path_raw);
raw_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mda',path_raw);
raw_subset_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17_subset.mda',path_raw);

if 1
%if (~exist(raw_mda_fname,'file'))
    fprintf('Loading raw data...\n');
    L=load(raw_mat_fname);
    raw=L.dl12_20151208_NNF_r1_tet16_17.channelData';
    fprintf('Writing raw data...\n');
    writemda(raw,raw_mda_fname);
    %spikespy({raw});
    raw=raw([1,3:6],(1e6+1):26e6);
    %raw=raw(7:10,(1e6+1):26e6);
    fprintf('Writing raw subset data...\n');
    writemda(raw,raw_subset_mda_fname);
end;

writemda(get_geometry,sprintf('%s/locations.mda',path0));

mscmd_bandpass_filter(raw_subset_mda_fname,[path0,'/pre1.mda'],opts_pre.o_filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],opts_pre.o_whiten);
mscmd_bandpass_filter([path0,'/pre2.mda'],[path0,'/pre3.mda'],opts_pre.o_filter2);

mscmd_detect([path0,'/pre3.mda'],[path0,'/detect0.mda'],opts.o_detect);

fprintf('Reading pre3...\n');
pre3=readmda([path0,'/pre3.mda']);

fprintf('Reading detect...\n');
detect0=readmda([path0,'/detect0.mda']);
times0=detect0(2,:);

fprintf('Extracting clips...\n');
clips0=ms_extract_clips(pre3,times0,opts.o_features.clip_size);
[M,T,NC]=size(clips0);

clusters=do_sorting(clips0,opts);

fprintf('Assembling templates...\n');
templates=zeros(M,T,0);
for j=1:length(clusters)
    templates=cat(3,templates,clusters{j}.template);
end;
fprintf('Writing...\n');
writemda(templates,[path0,'/templates.mda']);

fprintf('Assembling clips and clusters...\n');
clips=zeros(M,T,0);
clips_index=[];
all_inds=[];
labels=[];
ii=0;
for j=1:length(clusters)
    clips_index=[clips_index,ii]; ii=ii+length(clusters{j}.inds);
    clips=cat(3,clips,clips0(:,:,clusters{j}.inds));
    all_inds=[all_inds,clusters{j}.inds];
    labels=[labels,ones(1,length(clusters{j}.inds))*j];
end;

detect=detect0(:,all_inds);
times=detect(2,:);
[~,times_sort_inds]=sort(times);
detect=detect(:,times_sort_inds);
labels=labels(times_sort_inds);

clusters=zeros(3,length(all_inds));
clusters(1:2,:)=detect(1:2,:);
clusters(3,:)=labels;

fprintf('Writing...\n');
writemda(clips,[path0,'/clips.mda']);
writemda(clips_index,[path0,'/clips_index.mda']);
writemda(detect,[path0,'/detect.mda']);
writemda(clusters,[path0,'/clusters.mda']);

mscmd_cross_correlograms([path0,'/clusters.mda'],[path0,'/cross_correlograms.mda'],opts.o_cross_correlograms.max_dt);

mscmd_extract_clips([path0,'/pre1.mda'],[path0,'/clusters.mda'],[path0,'/clips_pre1.mda'],[path0,'/clips_index_pre1.mda'],opts.o_extract_clips);
mscmd_templates([path0,'/pre1.mda'],[path0,'/clusters.mda'],[path0,'/templates_pre1.mda'],opts.o_templates);

mscmd_templates(raw_subset_mda_fname,[path0,'/clusters.mda'],[path0,'/templates_raw.mda'],opts.o_templates);



end

function clusters=do_sorting(clips,opts)

[N,T,NC]=size(clips);
Tcenter=ceil((T+1)/2);

clusters={};

isosplit_opts.isocut_threshold=opts.o_cluster.ks_threshold;
isosplit_opts.K=opts.o_cluster.K_init;

fprintf('features...\n');
features=ms_event_features(clips,opts.o_features.num_features);
fprintf('isosplit...\n');
labels=isosplit(features,isosplit_opts);

K=max(labels);

fprintf('%d clusters found (size %d).\n',K,length(labels));

aa=opts.o_cluster.min_cluster_split_size;

for k=1:K
    indices_k=find(labels==k);
    clips0=clips(:,:,indices_k);
    if (size(clips0,3)>aa*2)
        peaks=squeeze(max(abs(clips0(:,Tcenter,:)),[],1));
        [~,sort_inds]=sort(peaks);
        inds_to_use=sort_inds(aa:end);
        clips1=clips0(:,:,inds_to_use);
        clusters2=do_sorting(clips1,opts);
        if (length(clusters2)>1)
            for j=1:length(clusters2)
                CC=clusters2{j};
                CC.inds=indices_k(inds_to_use(CC.inds));
                clusters{end+1}=CC;
            end;
        else
            CC.template=mean(clips1,3);
            CC.inds=indices_k;
            clusters{end+1}=CC;
        end;
    else
        CC.template=mean(clips0,3);
        CC.inds=indices_k;
        clusters{end+1}=CC;
    end;
end;

end

function opts_pre=get_preprocessing_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial bandpass filter options
opts_pre.o_filter.samplefreq=30000;
opts_pre.o_filter.freq_min=600;
opts_pre.o_filter.freq_max=4000;
opts_pre.o_filter.outlier_threshold=500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Whitening options
opts_pre.o_whiten=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-whitening bandpass filter options
opts_pre.o_filter2.samplefreq=30000;
opts_pre.o_filter2.freq_min=600;
opts_pre.o_filter2.freq_max=4000;
opts_pre.o_filter2.outlier_threshold=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function opts=get_sorting_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detection options
opts.o_detect.inner_window_width=15;
opts.o_detect.outer_window_width=1000;
opts.o_detect.threshold=2.5;
opts.o_detect.use_pre1=0;
opts.o_detect.individual_channels=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature extraction options
opts.o_features.num_features=4;
opts.o_features.clip_size=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clustering options
opts.o_cluster.ks_threshold=1.2;
opts.o_cluster.K_init=30;
opts.o_cluster.min_cluster_split_size=600;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster splitting options
% opts.o_split_clusters.num_features=3; 
% opts.o_split_clusters.clip_size=opts.o_features.clip_size;
% opts.o_split_clusters.use_pre1=1;
% opts.o_split_clusters.ks_threshold=opts.o_cluster.ks_threshold;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Template extraction options
opts.o_templates.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster consolidation options
%opts.o_consolidate.coincidence_threshold=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit options
opts.o_fit=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross-correlogram options
opts.o_cross_correlograms.max_dt=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clip extraction options
opts.o_extract_clips.clip_size=opts.o_templates.clip_size;
end

function L=get_geometry
L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];

end
