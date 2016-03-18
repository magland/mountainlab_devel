function firings_out=ms_outlier_scores_v1(X,firings,opts)
%MS_OUTLIER_SCORES - Compute outlier scores and store them in the 5th
%row of the firings array.
%
%Consider using mscmd_outlier_scores_v1
%
% Syntax:  [firings_out] = ms_outlier_scores_v1(X,firings,opts)
%
% Inputs:
%    X - MxN array of preprocessed raw data
%    firings - RxL array of firings containing times/labels etc (see docs)
%    opts.clip_size - the clip size in timepoints, e.g. 100
%    opts.shell_increment - Controls the definition of subshells, e.g. 0.5
%    opts.min_shell_size - Controls the definition of subshells, e.g. 100
%
% Outputs:
%    firings_out - modified version of firings with outlier scores stored
%    in the 5th row
%
% Other m-files required: ms_event_features, ms_extract_clips, ms_geometric_median
%
% See also: mscmd_outlier_scores_v1

% Author: Jeremy Magland
% Feb 2016; Last revision: 29-Feb-2016 (leap day)

if (nargin<1) test_outlier_scores_v1; return; end;

if (nargin<3) opts=struct; end;
if (~isfield(opts,'clip_size')) opts.clip_size=100; end;
if (~isfield(opts,'shell_increment')) opts.shell_increment=0.5; end;
if (~isfield(opts,'min_shell_size')) opts.min_shell_size=100; end;

L=size(firings,2);
[M,N]=size(X);

times=firings(2,:);
labels=firings(3,:);
peaks=firings(4,:);
scores=zeros(1,L); %outlier scores
K=max(labels);

clips=ms_extract_clips(X,times,opts.clip_size);

%Define random clips
interval=ceil(N/5000);
ttt=1:interval:N;
random_clips=ms_extract_clips(X,ttt,opts.clip_size);

for k=1:K
    inds_k=find(labels==k);
    peaks_k=peaks(inds_k);
    shells_k=define_shells(peaks_k,opts);
    for s=1:length(shells_k)
        inds_ks=inds_k(shells_k{s});
        clips_ks=clips(:,:,inds_ks);
        scores_ks=compute_outlier_scores(clips_ks,random_clips,opts);
        scores(inds_ks)=scores_ks;
    end;
end

firings_out=firings;
firings_out(5,:)=scores;

end

function scores=compute_outlier_scores(clips,random_clips,opts)
[M,T,L]=size(clips);
num_random_clips=size(random_clips,3);
template0=mean(clips,3);
weights=get_template_weights(template0,6);
random_clips_weighted=random_clips.*repmat(weights,1,1,num_random_clips);
clips_weighted=clips.*repmat(weights,1,1,size(clips,3));
template_weighted=template0.*weights;
diffs_weighted=clips_weighted-repmat(template_weighted,1,1,L);

% num_features=6;
% [FF,subspace]=ms_event_features(clips_weighted,num_features);
% for j=1:num_features
%     factor=sqrt(var(FF(j,:)));
%     subspace(:,:,j)=subspace(:,:,j)/factor;
% end;
% FF_random_clips=zeros(num_features,num_random_clips);
% for j=1:num_features
%     FF_random_clips(j,:)=squeeze(sum(sum(repmat(subspace(:,:,j),1,1,num_random_clips).*random_clips_weighted,1),2));
% end;
% FF_diffs=zeros(num_features,L);
% for j=1:num_features
%     FF_diffs(j,:)=squeeze(sum(sum(repmat(subspace(:,:,j),1,1,L).*diffs_weighted,1),2));
% end;
% vals1=sum(FF_diffs.^2,1);
% vals2=sum(FF_random_clips.^2,1);

vals1=reshape(sum(sum(diffs_weighted.^2,1),2),1,L);
vals2=reshape(sum(sum(random_clips_weighted.^2,1),2),1,num_random_clips);
mu0=mean(vals2); sigma0=sqrt(var(vals2));
scores=(vals1-mu0)/sigma0;
end

function W=get_template_weights(template,num_pix)
[M,T]=size(template);

aa=ifftshift(-floor(T/2):-floor(T/2)+T-1);
sig=num_pix;
kernel=exp(-0.5*aa.^2/sig^2);

fhat=fft(abs(template),[],2);
fhat=fhat.*repmat(kernel,M,1);
W=real(ifft(fhat,[],2));
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
peak_upper=peak_lower+opts.shell_increment;
while (peak_lower<=max_peak)
    inds1=find((peak_lower<=peaks)&(peaks<peak_upper));
    ct1=length(inds1);
    ct2=length(find(peaks>=peak_upper));
    if (peak_upper>max_peak)
        shells{end+1}=inds1;
        peak_lower=peak_upper;
    elseif ((ct1>=opts.min_shell_size)&&(ct2>=opts.min_shell_size))
        shells{end+1}=inds1;
        peak_lower=peak_upper;
        peak_upper=peak_lower+opts.shell_increment;
    else
        peak_upper=peak_upper+opts.shell_increment;
    end;
end;

end

function test_outlier_scores_v1
mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../demo/demo_sort_002/output_tetrode1'];
X=readmda([path0,'/pre2.mda']);
firings=readmda([path0,'/firings.mda']);
firings2=ms_outlier_scores_v1(X,firings);
writemda(firings2,'tmp_firings2.mda');
mv.pre=[path0,'/pre2.mda'];
mv.firings='tmp_firings2.mda';
mv.sampling_freq=30000;
ms_mountainview(mv);
end

