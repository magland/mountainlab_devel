function example_validation

addpath('mountainview/src/spikespy/matlab');
addpath('processing');
addpath('msutils');

mfile_path=fileparts(mfilename('fullpath'));

path_raw=sprintf('%s/../example_data',mfile_path);
path0=sprintf('%s/output',mfile_path);
path1=sprintf('%s/output_reversed',mfile_path);

if (~exist(path0,'dir')) mkdir(path0); end;
if (~exist(path1,'dir')) mkdir(path1); end;

opts=get_sorting_options;

if (~exist([path_raw,'/ms11d45.dat'],'file'))
    fprintf('File not found: %s\n',[path0,'/ms11d45.dat']);
    return;
end;

if ((~exist([path0,'/adjacency.mda'],'file'))||(~exist([path0,'/locations.mda'],'file')))
    locations=get_frank_lab_locations;
    AM=ms_adjacency_matrix(locations,2);
    writemda(AM,[path0,'/adjacency.mda']);
    writemda(locations,[path0,'/locations.mda']);
end;

mscmd_extract([path_raw,'/ms11d45.dat'],[path0,'/raw.mda'],opts.o_extract);
mscmd_bandpass_filter([path0,'/raw.mda'],[path0,'/filt.mda'],opts.o_filter);
mscmd_whiten([path0,'/filt.mda'],[path0,'/filt_white.mda'],opts.o_whiten);
mscmd_bandpass_filter([path0,'/filt_white.mda'],[path0,'/filt2_white.mda'],opts.o_filter2);

do_sorting(path0,opts);

reverse_noise(path0,path1);

do_sorting(path1,opts);

fprintf('MAPPING/MATCHING CLUSTERS\n');
C1=readmda(sprintf('%s/clusters.mda',path0));
C2=readmda(sprintf('%s/clusters.mda',path1));
C2_mapped=map_clusters(C1,C2);
writemda(C2_mapped,sprintf('%s/clusters_mapped.mda',path1));

example_validation_view;

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
opts.o_filter.freq_min=600;
opts.o_filter.freq_max=4000;
opts.o_filter.outlier_threshold=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_whiten.ncomp=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_filter2.samplefreq=30000;
opts.o_filter2.freq_min=600;
opts.o_filter2.freq_max=4000;
opts.o_filter2.outlier_threshold=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_detect.inner_window_width=15;
opts.o_detect.outer_window_width=100000;
opts.o_detect.threshold=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_features.num_features=6;
opts.o_features.clip_size=60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_cluster=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_split_clusters.num_features=3; 
opts.o_split_clusters.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_templates.clip_size=opts.o_features.clip_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_consolidate.coincidence_threshold=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_cross_correlograms.max_dt=1500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.o_extract_clips.clip_size=opts.o_templates.clip_size;

end

function do_sorting(output_path,opts)

timerA=tic;

path0=output_path;

if ((~exist([path0,'/adjacency.mda'],'file'))||(~exist([path0,'/locations.mda'],'file')))
    locations=get_frank_lab_locations;
    AM=ms_adjacency_matrix(locations,2);
    writemda(AM,[path0,'/adjacency.mda']);
    writemda(locations,[path0,'/locations.mda']);
end;

mscmd_detect([path0,'/filt2_white.mda'],[path0,'/detect.mda'],opts.o_detect);
mscmd_features([path0,'/filt2_white.mda'],[path0,'/detect.mda'],[path0,'/adjacency.mda'],[path0,'/features.mda'],opts.o_features);
%remove_zero_cluster([path0,'/features.mda'],[path0,'/features2.mda']);
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
mscmd_cross_correlograms([path0,'/clusters.mda'],[path0,'/cross-correlograms.mda'],opts.o_cross_correlograms.max_dt);

fprintf('Elapsed time: %g sec\n',toc(timerA));

end

function C2=map_clusters(C1,C2)

[CM,map12,map21]=confusion_matrix(C1(2,:),C1(3,:),C2(2,:),C2(3,:));

used=zeros(1,length(map21));
used(map12(find(map12>0)))=1;
unused_inds=find(used==0);
new_order=[map12,unused_inds];
[~,new_order_inv]=sort(new_order);
C2(3,:)=new_order_inv(C2(3,:));

end

function reverse_noise(path0,path1)

fprintf('REVERSING NOISE...\n');

X=readmda([path0,'/filt2_white.mda']);
C=readmda([path0,'/clusters.mda']);
W=readmda([path0,'/templates0_filt2_white.mda']);

X2=X;
[M,T,K]=size(W);
[~,N]=size(X);

tt1=-floor(T/2); tt2=tt1+T-1;

L=size(C,2);
C_times=C(2,:);
C_labels=C(3,:);
for j=1:L
    if (mod(j,10000)==1) fprintf('%g%% ',floor(j/L*100)); end;
    t0=C_times(j);
    k=C_labels(j);
    X2(:,t0+tt1:t0+tt2)=X2(:,t0+tt1:t0+tt2)-2*W(:,:,k);
end;
X2=-X2;

writemda(X2,[path1,'/filt2_white.mda']);

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