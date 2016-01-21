function test_2016_01_21_validate

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

mfile_path=fileparts(mfilename('fullpath'));
addpath(sprintf('%s/..',mfile_path));
outputdir_path=sprintf('%s/output',mfile_path);
if (~exist(outputdir_path,'dir')) mkdir(outputdir_path); end;
datfile_path=sprintf('%s/../raw/ms11d45.dat',mfile_path);
path0=sprintf('%s/../experiment1/output',mfile_path);


locations=get_frank_lab_locations;
opts.adjacency=ms_adjacency_matrix(locations,2);
writemda(locations,[outputdir_path,'/locations.mda']);

if 0
    reverse_noise(sprintf('%s/pre3.mda',path0),sprintf('%s/clusters.mda',path0),sprintf('%s/templates.mda',path0),sprintf('%s/pre3.mda',outputdir_path));
end;

mountainsort_v1('',outputdir_path,opts);

clusters1=sprintf('%s/clusters.mda',path0);
clusters2=sprintf('%s/clusters.mda',outputdir_path);
mscmd_confusion_matrix(clusters1,clusters2,sprintf('%s/validation_matrix.mda',outputdir_path),3);

end

function reverse_noise(path0,clusters_path,templates_path,output_path)

fprintf('REVERSING NOISE...\n');

X=readmda(path0);
C=readmda(clusters_path);
W=readmda(templates_path);

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

writemda(X2,output_path);

end