function compare_ms_kk

mfile_path=fileparts(mfilename('fullpath'));

addpath([mfile_path,'/../processing']);
addpath([mfile_path,'/../msutils']);

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
end;
CC_kk=readmda(sprintf('%s/clusters.mda',output_kk_path));

% Find the mapping from kk to ms
fprintf('Compute cluster mapping from kk to ms...\n');
[mapping_kk,CM]=compute_cluster_mapping(sprintf('%s/clusters.mda',output_kk_path),sprintf('%s/clusters.mda',output_ms_path),sprintf('%s/confusion_matrix.mda',output_kk_path));
writemda(mapping_kk,sprintf('%s/mapping_kk.mda',output_kk_path));
writemda(CM,sprintf('%s/CM.mda',output_kk_path));
CC_kk(3,:)=mapping_kk(CC_kk(3,:));
inds=find(CC_kk(3,:)<=K_ms);
CC_kk=CC_kk(:,inds);
writemda(CC_kk,sprintf('%s/clusters_mapped.mda',output_kk_path));

figure; imagesc(CM');
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

fprintf('Computing cross-correlograms...\n');
mscmd_cross_correlograms([output_kk_path,'/clusters.mda'],[output_kk_path,'/cross-correlograms.mda'],opts.o_cross_correlograms.max_dt);

fprintf('Computing cross-correlograms-mapped...\n');
mscmd_cross_correlograms([output_kk_path,'/clusters_mapped.mda'],[output_kk_path,'/cross-correlograms-mapped.mda'],opts.o_cross_correlograms.max_dt);

view_ms(output_ms_path,output_kk_path);
%view_kk(output_ms_path,output_kk_path);
%view_compare_labels(output_ms_path,output_kk_path);

end

function view_ms(output_ms_path,output_kk_path)

mfile_path=fileparts(mfilename('fullpath'));
exe_fname=sprintf('%s/../mountainview/bin/mountainview',mfile_path);

cmd='';
cmd=[cmd,sprintf('%s ',exe_fname)];
cmd=[cmd,sprintf('--raw=%s/filt2_white.mda ',output_ms_path)];
cmd=[cmd,sprintf('--templates=%s/templates0_filt2_white.mda ',output_ms_path)];
cmd=[cmd,sprintf('--clips=%s/clips_filt2_white.mda ',output_ms_path)];
cmd=[cmd,sprintf('--primary-channels=%s/load_channels0.mda ',output_ms_path)];

cmd=[cmd,sprintf('--cluster=%s/clusters.mda ',output_ms_path)];
cmd=[cmd,sprintf('--locations=%s/locations.mda ',output_ms_path)];
cmd=[cmd,sprintf('--cross-correlograms=%s/cross-correlograms.mda ',output_ms_path)];

cmd=[cmd,sprintf('--clips-index=%s/clips_filt_index.mda',output_ms_path)];

fprintf('%s\n',cmd);
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
cmd=[cmd,sprintf('--cross-correlograms=%s/cross-correlograms-mapped.mda ',output_kk_path)];

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

function [mapping_kk,CM2]=compute_cluster_mapping(clusters_kk_path,clusters_ms_path,confusion_matrix_path)

mscmd_confusion_matrix(clusters_kk_path,clusters_ms_path,confusion_matrix_path,3);
CM=readmda(confusion_matrix_path);
CM=CM(1:end-1,1:end-1);
[K1,K2]=size(CM);

CM2=CM;
for j1=1:K1
    for j2=1:K2
        if ((sum(CM(:,j2))+sum(CM(j1,:)))~=0)
            CM2(j1,j2)=2*CM(j1,j2)/(sum(CM(:,j2))+sum(CM(j1,:)));
        end;
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
    out(inds_unused)=(size(CM,2)+1):(size(CM,2)+length(inds_unused));
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

function remove_zero_cluster(features_in_path,features_out_path)

fprintf('*** REMOVE ZERO CLUSTER ***\n');
tmp_in=readmda(features_in_path);
channels=tmp_in(1,:);
times=tmp_in(2,:);
features=tmp_in(3:end,:);
Nf=size(features,1);

M=max(channels);
channels_out=[];
times_out=[];
features_out=zeros(Nf,0);
for ch=1:M
    fprintf('Channel %d\n');
    inds0=find(channels==ch);
    times0=times(inds0);
    features0=features(:,inds0);
    okay_inds=find_nonzero_vectors(features0,25,1);
    fprintf('Using %d/%d events\n',length(okay_inds),length(inds0));
    channels_out=[channels_out,ones(1,length(okay_inds))*ch];
    times_out=[times_out,times0(okay_inds)];
    features_out=cat(2,features_out,features0(:,okay_inds));
end;

fprintf('Using a total of %d/%d events\n',length(times_out),length(times));

tmp_out=zeros(2+Nf,length(times_out));
tmp_out(1,:)=channels_out;
tmp_out(2,:)=times_out;
tmp_out(3:end,:)=features_out;
writemda(tmp_out,features_out_path);

end

function inds=find_nonzero_vectors(V,num_neighbor_pts,verbose)

cutoff=3;

[M,N]=size(V);
sigma_guess=sqrt(var(V(:)));
norms=sqrt(sum(V.^2,1));

fprintf('knnsearch...\n');
[ids,ds]=knnsearch(V',V','K',num_neighbor_pts);

% Volume of sphere is (pi)^[M/2] / gamma((M+2)/2) * R^M
% where R=ds(:,end)
log_est_densities=log(num_neighbor_pts)-(M/2)*log(pi)+gammaln((M+2)/2)-M*log(ds(:,end)');

fprintf('estimate_sigma...\n');
[est_sigma,est_log_density_at_zero]=estimate_sigma(M,norms,log_est_densities,sigma_guess);

% PDF of gaussian: (2*pi)^(-M/2)*sigma^M*exp(-norm^2/(2*sigma^2))
%log_expected_densities=log(N)-M/2*log(2*pi)-M*log(est_sigma)-1/2*norms.^2/est_sigma^2;
log_expected_densities=est_log_density_at_zero-1/2*(norms).^2/est_sigma^2;

if (verbose)
figure; plot(log_expected_densities,log_est_densities,'b.');
xlabel('Log expected densities');
ylabel('Log estimated densities');
hold on;
plot(xlim,xlim,'r');
end;

scores=log_est_densities-log_expected_densities;

if (verbose)
figure;
plot(norms,scores,'b.');
hold on;
plot(norms(scores>cutoff),scores(scores>cutoff),'r.');
xlabel('Norm');
ylabel('Score');
end;

if (verbose)
figure; set(gcf,'position',[100,1200,1000,500]);
subplot(1,2,1);
bins=linspace(min(scores),max(scores),1000);
hist(scores,bins); hold on;
hist(scores(scores>cutoff),bins); hold on;
xlabel('Scores');

subplot(1,2,2);
bins=linspace(min(norms),max(norms),1000);
hist(norms,bins); hold on;
hist(norms(scores>cutoff),bins); hold on;
xlabel('Norms');
end;

inds=find(scores>cutoff);

end

function [est_sigma,est_log_density_at_zero]=estimate_sigma(M,norms,log_est_densities,sigma_guess)
%figure; plot(norms.^2,log_est_densities,'b.');
%xlabel('norms^2');
%ylabel('log_est_densities');
inds=find(norms<=sigma_guess*sqrt(M));
[slope,intercept]=estimate_slope_intercept(norms(inds).^2,log_est_densities(inds));
est_sigma=sqrt(-1/(2*slope));
est_log_density_at_zero=intercept;
end

function [slope,intercept]=estimate_slope_intercept(X,Y)
coeffs=polyfit(X,Y,1);
slope=coeffs(1);
intercept=coeffs(2);
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
