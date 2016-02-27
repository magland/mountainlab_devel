function [firings_path,pre_path]=sort_003_multichannel(raw_path,output_path,sort_opts)
%SORT_003_MULTICHANNEL - Version 003 of sorting based on shell method and isosplit2
%
% Syntax:  firings_path=sort_003_multichannel(raw_path,output_path,sort_opts)
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

if nargin<1 test_sort_003_multichannel; return; end;

def_sort_opts.clip_size=120;
def_sort_opts.filter.samplefreq=30000;
def_sort_opts.filter.freq_min=300;
def_sort_opts.filter.freq_max=6000;
def_sort_opts.filter.outlier_threshold=500;
def_sort_opts.detect.detect_threshold=3;
def_sort_opts.detect.detect_interval=15;
def_sort_opts.detect.individual_channels=1;
def_sort_opts.plausibility_threshold=3;
def_sort_opts.detectability_threshold=9;
%def_sort_opts.min_cluster_size=10;
%def_sort_opts.detect.individual_channels=1;
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

%%%% Clustering
mscmd_branch_cluster_v1([path0,'/pre2.mda'],[path0,'/detect.mda'],'',[path0,'/firings1.mda']);

%%%% Pruning
ooo.detectability_threshold=sort_opts.detectability_threshold;
ooo.clip_size=sort_opts.clip_size;
mscmd_remove_duplicates([path0,'/firings1.mda'],[path0,'/firings2.mda']);
%remove_noisy_subclusters([path0,'/pre2.mda'],[path0,'/firings2.mda'],[path0,'/firings3.mda'],ooo);

%%%% Copying
mscmd_copy([path0,'/firings2.mda'],[path0,'/firings.mda']);
%mscmd_copy([path0,'/firings3.mda'],[path0,'/firings.mda']);

firings_path=[path0,'/firings.mda'];
pre_path=[path0,'/pre2.mda'];

end

function score=compute_detectability_score(template,rand_clips,ch)
[M,T,NC]=size(rand_clips);
%subtract the mean clip
mean_clip=mean(rand_clips,3);

rand_clips=rand_clips-repmat(mean(rand_clips,2),1,T,1);
template=template-repmat(mean(template,2),1,T);
weights=get_template_weights(template);
Wtemplate=template.*weights;
Wrand_clips=rand_clips.*repmat(weights,1,1,NC);
Vtemplate=reshape(Wtemplate,M*T,1);
Vrand_clips=reshape(Wrand_clips,M*T,NC);
covmat=V
end

function score=compute_detectability_score_old(template,rand_clips)
[M,T,L]=size(rand_clips);
inner_products=squeeze(sum(sum(rand_clips.*repmat(template,1,1,L),1),2));
inner_product_0=sum(sum(template.*template,1),2);
mu=mean(inner_products);
sigma=sqrt(var(inner_products));
score=(inner_product_0-mu)/sigma;
%figure; hist(inner_products,1000); hold on; plot([inner_product_0,inner_product_0],ylim,'linewidth',6);
end

function remove_noisy_subclusters(pre_path,firings_path,firings_out_path,opts)
firings=readmda(firings_path);
L=size(firings,2);
pre=readmda(pre_path);
channels=firings(1,:); times=firings(2,:); labels=firings(3,:); peaks=firings(4,:);
K=max(labels);
inds_to_remove=[];
num_subclusters_removed=0;
for k=1:K
    inds=find(labels==k);
    channel=channels(inds(1));
    random_clips=extract_random_clips(pre,channel,0,5000,opts.clip_size); %important to have low threshold!
    shells=define_shells(peaks(inds),struct('increment',0.5,'min_per_shell',100));
    for s=1:length(shells)
        inds0=inds(shells{s});
        clips0=ms_extract_clips(pre,times(inds0),opts.clip_size);
        template0=compute_geometric_median_template(clips0);
        score0=compute_detectability_score(template0,random_clips)
        if (score0<opts.detectability_threshold)
            inds_to_remove=[inds_to_remove,inds0];
            num_subclusters_removed=num_subclusters_removed+1;
        end;
    end;
end;
fprintf('Removing %d events in %d noisy subclusters...\n',length(inds_to_remove),num_subclusters_removed);
to_use=ones(1,L);
to_use(inds_to_remove)=0;
firings=firings(:,find(to_use));

%remap the labels
labels=firings(3,:);
label_map=zeros(1,K);
kk=1;
for k=1:K
    if (length(find(labels==k))>0)
        label_map(k)=kk;
        kk=kk+1;
    end;
end;
labels(find(labels~=0))=label_map(labels(find(labels~=0)));
firings(3,:)=labels;

fprintf('Using %d of %d clusters...\n',max(labels),K);

writemda(firings,firings_out_path);

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

function clips=extract_random_clips(pre,channel,threshold,num,clip_size)
[M,N]=size(pre);
%times=ms_detect(pre,struct('detect_threshold',threshold,'detect_interval',clip_size,'clip_size',clip_size));
times=clip_size+randsample(N-2*clip_size,num,true)';
if (length(times)>num)
    incr=floor(length(times)/num);
    times=times((1:num)*incr);
end;
clips=ms_extract_clips(pre,times,clip_size);
end

function shells=define_shells(peaks,opts)
shells={};

%negatives and positives
if (min(peaks)<0)
    inds_neg=find(peaks<0);
    inds_pos=find(peaks>=0);
    shells_neg=define_shells(-peaks(inds_neg),opts);
    shells_pos=define_shells(peaks(inds_pos),opts);
    for ii=length(shells_neg):-1:1
        shells{end+1}=inds_neg(shells_neg{ii});
    end;
    for ii=1:length(shells_pos)
        shells{end+1}=inds_pos(shells_pos{ii});
    end;
    return;
end;

%only positives
max_peak=max(peaks);
peak_lower=0;
peak_upper=peak_lower+opts.increment;
while (peak_lower<=max_peak)
    inds1=find((peak_lower<=peaks)&(peaks<peak_upper));
    ct1=length(inds1);
    ct2=length(find(peaks>=peak_upper));
    if (peak_upper>max_peak)
        shells{end+1}=inds1;
        peak_lower=peak_upper;
    elseif ((ct1>=opts.min_per_shell)&&(ct2>=opts.min_per_shell))
        shells{end+1}=inds1;
        peak_lower=peak_upper;
        peak_upper=peak_lower+opts.increment;
    else
        peak_upper=peak_upper+opts.increment;
    end;
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


function test_sort_003_multichannel

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
[firings_path,pre_path]=sort_003_multichannel([path0,'/pre0.mda'],path0);

%%%% View output
mv.mode='overview2';
mv.raw=pre_path;
mv.firings=firings_path'
ms_mountainview(mv);

end


