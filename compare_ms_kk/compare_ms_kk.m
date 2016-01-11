function compare_ms_kk

mfile_path=fileparts(mfilename('fullpath'));
raw_path=sprintf('%s/../example_data/ms11d45.dat',mfile_path);
kwik_path=sprintf('%s/../example_data/ms11d45.kwik',mfile_path);
output_ms_path=sprintf('%s/output_ms',mfile_path);
output_kk_path=sprintf('%s/output_kk',mfile_path);
if ~exist(output_ms_path,'dir') mkdir(output_ms_path); end;
if ~exist(output_kk_path,'dir') mkdir(output_kk_path); end;
preprocessed_path=sprintf('%s/filt2_white.mda',output_ms_path);

% Do the MountainSort sorting
fprintf('MountainSort sorting...\n');
opts=get_sorting_options;
do_sorting(raw_path,output_ms_path,opts);
CC_ms=readmda(sprintf('%s/clusters.mda',output_ms_path));
K_ms=max(CC_ms(3,:));

% Read the KlustaKwik data
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

% Find the mapping from kk to ms
if 0
    fprintf('Compute cluster mapping from kk to ms...\n');
    [mapping_kk,CM]=compute_cluster_mapping(CC_kk(2,:),CC_kk(3,:),CC_ms(2,:),CC_ms(3,:));
    writemda(mapping_kk,sprintf('%s/mapping_kk.mda',output_kk_path));
    writemda(CM,sprintf('%s/CM.mda',output_kk_path));
    CC_kk(3,:)=mapping_kk(CC_kk(3,:));
    inds=find(CC_kk(3,:)<=K_ms);
    CC_kk=CC_kk(:,inds);
    writemda(CC_kk,sprintf('%s/clusters_mapped.mda',output_kk_path));
    figure; imagesc(CM');
end;
selected_clusters=[8, 25, 128, 19, 56, 100, 83, 98, 99, 65, 103, 94, 107, 82, 93, 95, 80, 104];
mapping_kk=readmda(sprintf('%s/mapping_kk.mda',output_kk_path));
disp('mapping kk->ms:');
ct=0;
for j=1:K_ms
    ind=find(mapping_kk==j);
    if (length(ind)>0)
        if (length(find(selected_clusters==ind))>0)
            fprintf('[%d->%d]   ',ind,j);
        else
            fprintf(' %d->%d    ',ind,j);
        end;
        ct=ct+1;
        if (mod(ct,6)==0) fprintf('\n'); end;
    end;
end;
fprintf('\n');

CC_ms=readmda(sprintf('%s/clusters.mda',output_ms_path));
to_keep=[];
for j=1:length(selected_clusters);
    to_keep=[to_keep,find(CC_ms(3,:)==mapping_kk(selected_clusters(j)))];
end;
CC_ms=CC_ms(:,to_keep);
writemda(CC_ms,sprintf('%s/clusters_ms_subset.mda',output_kk_path));
CC_kk=readmda(sprintf('%s/clusters.mda',output_kk_path));
to_keep=[];
for j=1:length(selected_clusters);
    to_keep=[to_keep,find(CC_kk(3,:)==selected_clusters(j))];
end;
CC_kk=CC_kk(:,to_keep);
writemda(CC_kk,sprintf('%s/clusters_kk_subset.mda',output_kk_path));

mscmd_templates(preprocessed_path,sprintf('%s/clusters_mapped.mda',output_kk_path),sprintf('%s/templates_mapped.mda',output_kk_path),opts.o_templates);

%fprintf('Computing cross-correlograms...\n');
%[cc_kk,CCmda_kk]=ms_cross_correlograms(CC_kk(2,:),CC_kk(3,:),opts.o_cross_correlograms.max_dt);
%writemda(CCmda,sprintf('%s/cross-correlograms.mda',output_kk_path));

view_ms(output_ms_path,output_kk_path);
view_kk(output_ms_path,output_kk_path);
view_compare_labels(output_ms_path,output_kk_path);

end

function view_ms(output_ms_path,output_kk_path)

mfile_path=fileparts(mfilename('fullpath'));
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s ',exe_fname)];
cmd=[cmd,sprintf('--raw=%s/filt2_white.mda ',output_ms_path)];
cmd=[cmd,sprintf('--templates=%s/templates0_filt2_white.mda ',output_ms_path)];
cmd=[cmd,sprintf('--clips=%s/clips_filt2_white.mda ',output_ms_path)];

cmd=[cmd,sprintf('--cluster=%s/clusters.mda ',output_ms_path)];
cmd=[cmd,sprintf('--locations=%s/locations.mda ',output_ms_path)];
cmd=[cmd,sprintf('--cross-correlograms=%s/cross-correlograms.mda ',output_ms_path)];

system(sprintf('%s &',cmd));

end

function view_kk(output_ms_path,output_kk_path)

mfile_path=fileparts(mfilename('fullpath'));
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s ',exe_fname)];
cmd=[cmd,sprintf('--raw=%s/filt2_white.mda ',output_ms_path)];
cmd=[cmd,sprintf('--templates=%s/templates_mapped.mda ',output_kk_path)];

cmd=[cmd,sprintf('--cluster=%s/clusters_mapped.mda ',output_kk_path)];
cmd=[cmd,sprintf('--locations=%s/locations.mda ',output_kk_path)];
%cmd=[cmd,sprintf('--cross-correlograms=%s/cross-correlograms.mda ',output_kk_path)];

system(sprintf('%s &',cmd));

end

function view_compare_labels(output_ms_path,output_kk_path)

mfile_path=fileparts(mfilename('fullpath'));
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s --mode=compare_labels ',exe_fname)];

cmd=[cmd,sprintf('--raw=%s/filt2_white.mda ',output_ms_path)];
cmd=[cmd,sprintf('--cluster=%s/clusters.mda ',output_ms_path)];
cmd=[cmd,sprintf('--cluster2=%s/clusters_kk_subset.mda ',output_kk_path)];

fprintf('%s\n',cmd);
system(sprintf('%s &',cmd));

end

function [mapping_kk,CM2]=compute_cluster_mapping(times_kk,labels_kk,times_ms,labels_ms)

[CM,map12,map21]=confusion_matrix(times_kk,labels_kk,times_ms,labels_ms);
CM=CM(1:end-1,1:end-1);
[K1,K2]=size(CM);

CM2=CM;
for j1=1:K1
    for j2=1:K2
        CM2(j1,j2)=2*CM(j1,j2)/(sum(CM(:,j2))+sum(CM(j1,:)));
    end;
end;

mapping_kk=get_optimal_mapping(CM2);

end

function out=get_optimal_mapping(CM)

inds1=[];
inds2=[];
while (max(CM(:)>0))
    [~,ind0]=max(CM(:));
    [i1,i2]=ind2sub(size(CM),ind0);
    inds1=[inds1,i1];
    inds2=[inds2,i2];
    CM(i1,:)=0;
    CM(:,i2)=0;
end;

out=zeros(1,size(CM,1));
out(inds1)=inds2;
inds_unused=find(out==0);
if (length(inds_unused)>0)
    out(inds_unused)=size(CM,2)+1:size(CM,1);
end;

end

function opts=get_sorting_options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_extract.num_channels=72;
opts.o_extract.channels=[37:52,68,69];
%opts.o_extract.t1=0; o_extract.t2=1e6;
opts.o_extract.t1=0; opts.o_extract.t2=19e6;
%opts.o_extract.t1=0; o_extract.t2=5e7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_filter.samplefreq=30000;
opts.o_filter.freq_min=300;
opts.o_filter.freq_max=4000;
opts.o_filter.outlier_threshold=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_whiten.ncomp=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_filter2.samplefreq=30000;
opts.o_filter2.freq_min=600;
opts.o_filter2.freq_max=10000;
opts.o_filter2.outlier_threshold=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_detect.inner_window_width=40;
opts.o_detect.outer_window_width=100000;
opts.o_detect.threshold=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_features.num_features=6;
opts.o_features.clip_size=200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_cluster=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_split_clusters.num_features=3; 
opts.o_split_clusters.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_templates.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_consolidate.compare_threshold=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_cross_correlograms.max_dt=1500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_extract_clips.clip_size=opts.o_templates.clip_size;

end

function do_sorting(raw_path,output_path,opts)

timerA=tic;

path0=output_path;

if ((~exist([path0,'/adjacency.mda'],'file'))||(~exist([path0,'/locations.mda'],'file')))
    locations=get_frank_lab_locations;
    AM=ms_adjacency_matrix(locations,2);
    writemda(AM,[path0,'/adjacency.mda']);
    writemda(locations,[path0,'/locations.mda']);
end;

mscmd_extract(raw_path,[path0,'/raw.mda'],opts.o_extract);
mscmd_bandpass_filter([path0,'/raw.mda'],[path0,'/filt.mda'],opts.o_filter);
mscmd_whiten([path0,'/filt.mda'],[path0,'/filt_white.mda'],opts.o_whiten);
mscmd_bandpass_filter([path0,'/filt_white.mda'],[path0,'/filt2_white.mda'],opts.o_filter2);

mscmd_detect([path0,'/filt2_white.mda'],[path0,'/detect.mda'],opts.o_detect);
mscmd_features([path0,'/filt2_white.mda'],[path0,'/detect.mda'],[path0,'/adjacency.mda'],[path0,'/features.mda'],opts.o_features);
mscmd_cluster([path0,'/features.mda'],[path0,'/cluster.mda'],opts.o_cluster);

mscmd_split_clusters([path0,'/filt2_white.mda'],[path0,'/cluster.mda'],[path0,'/cluster2.mda'],opts.o_split_clusters);
mscmd_templates([path0,'/filt2_white.mda'],[path0,'/cluster2.mda'],[path0,'/templates.mda'],opts.o_templates);
mscmd_consolidate([path0,'/cluster2.mda'],[path0,'/templates.mda'],[path0,'/cluster0.mda'],[path0,'/templates0.mda'],[path0,'/load_channels0.mda'],opts.o_consolidate);

mscmd_fit([path0,'/filt2_white.mda'],[path0,'/cluster0.mda'],[path0,'/templates0.mda'],[path0,'/clusters.mda']);
%writemda(readmda([path0,'/cluster0.mda']),[path0,'/clusters.mda']);

mscmd_templates([path0,'/raw.mda'],[path0,'/clusters.mda'],[path0,'/templates0_raw.mda'],opts.o_templates);
mscmd_templates([path0,'/filt2_white.mda'],[path0,'/clusters.mda'],[path0,'/templates0_filt2_white.mda'],opts.o_templates);

if (exist([path0,'/filt.mda']))
    mscmd_extract_clips([path0,'/filt.mda'],[path0,'/clusters.mda'],[path0,'/clips_filt.mda'],[path0,'/clips_filt_index.mda'],opts.o_extract_clips);
end;
mscmd_extract_clips([path0,'/filt2_white.mda'],[path0,'/clusters.mda'],[path0,'/clips_filt2_white.mda'],[path0,'/clips_filt2_white_index.mda'],opts.o_extract_clips);

fprintf('Computing cross-correlograms...\n');
clusters=readmda([path0,'/clusters.mda']);
[CC,CCmda]=ms_cross_correlograms(clusters(2,:),clusters(3,:),opts.o_cross_correlograms.max_dt);
writemda(CCmda,[path0,'/cross-correlograms.mda']);

fprintf('Elapsed time: %g sec\n',toc(timerA));

end

function L=get_frank_lab_locations

L=[...
0.0,0.0;...
-0.5,1.0;...
0.5,1.0;...
-1.0,2.0;...
1.0,2.0;...
-1.0,3.0;...
1.0,3.0;...
-1.0,4.0;...
1.0,4.0;...
-1.0,5.0;...
1.0,5.0;...
-1.0,6.0;...
1.0,6.0;...
-1.0,7.0;...
1.0,7.0;...
-1.0,8.0;...
1.0,8.0;...
-1.0,9.0;...
];

end
