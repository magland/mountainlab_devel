function firings_out=ms_remove_noise_subclusters(pre,firings,opts)
%MS_REMOVE_NOISE_SUBCLUSTERS - Remove noise subclusters by splitting clusters
%into shells and comparing the subcluster template waveforms with the
%expected waveform shape in the situation of pure noise detected via
%amplitude threshold.
%
%Better to use mscmd_remove_noise_subclusters
%
% Syntax:  [firings_out] = ms_remove_noise_subclusters(pre,firings,opts)
%
% Inputs:
%    pre - MxN array of preprocessed raw data
%    firings - RxL array of firings containing times/labels etc (see docs)
%    opts.clip_size - the clip size in timepoints, e.g. 100
%    opts.detectability_threshold - the threshold for accepting a
%             subcluster as non-noise. For example use 4.
%    opts.shell_increment - Controls the definition of subshells, e.g. 0.5
%    opts.min_shell_size - Controls the definition of subshells, e.g. 100
%
% Outputs:
%    firings2 - modified version of firings with events in noise
%             subclusters remove and labels re-mapped based on entire clusters that
%             have been removed.
%
% Other m-files required: ms_event_features, ms_extract_clips, ms_geometric_median
%
% See also: mscmd_remove_noise_subclusters

% Author: Jeremy Magland
% Feb 2016; Last revision: 29-Feb-2016 (leap day)
L=size(firings,2);
T=opts.clip_size;
Tmid=floor((T+1)/2);
channels=firings(1,:); times=firings(2,:); labels=firings(3,:); peaks=firings(4,:);
clips=ms_extract_clips(pre,times,opts.clip_size);
K=max(labels);
subclusters={};
for k=1:K
    fprintf('.');
    inds_k=find(labels==k);
    if (length(inds_k)>0)
        channel=channels(inds_k(1));
        subclusters0=compute_subcluster_detectability_scores(pre,clips(:,:,inds_k),channel,opts);
        for ii=1:length(subclusters0)
            subcluster=subclusters0{ii};
            subcluster.inds=inds_k(subcluster.inds);
            subcluster.label=k;
            subclusters{end+1}=subcluster;
        end;
    end;
end;
fprintf('\n');
num_subclusters=length(subclusters);
subcluster_peaks=zeros(1,num_subclusters);
subcluster_scores=zeros(1,num_subclusters);
for ii=1:num_subclusters
    subcluster_peaks(ii)=subclusters{ii}.peak;
    subcluster_scores(ii)=subclusters{ii}.detectability_score;
end;
coeffs=polyfit(abs(subcluster_peaks),subcluster_scores,1);
subcluster_scores=subcluster_scores/coeffs(1);
figure; plot(abs(subcluster_peaks),subcluster_scores,'b.','MarkerSize',8);
xlabel('Abs. subcluster peak'); ylabel('Subcluster detectability score');
labels_to_use=zeros(1,K);
events_to_use=zeros(1,L);
for ii=1:num_subclusters
    if (subcluster_scores(ii)>=opts.detectability_threshold)
        labels_to_use(subclusters{ii}.label)=1;
        events_to_use(subclusters{ii}.inds)=1;
    end;
end;
labels_map=zeros(1,K);
labels_map(find(labels_to_use))=1:length(find(labels_to_use));
firings_out=firings(:,find(events_to_use));
labels=firings_out(3,:);
labels(labels~=0)=labels_map(labels(labels~=0));
firings_out(3,:)=labels;
fprintf('Using %d/%d events in %d/%d clusters\n',size(firings_out,2),L,max(labels),K);

end


function subclusters=compute_subcluster_detectability_scores(pre,clips,channel,opts)
[M,T,L]=size(clips);
Tmid=floor((T+1)/2);
peaks=reshape(clips(channel,Tmid,:),1,L);
shells=define_shells(peaks,opts);
noise_shape=estimate_noise_shape(pre,T,channel);
subclusters={};
for s=1:length(shells)
    inds_s=shells{s};
    clips_s=clips(:,:,inds_s);
    subtemplate=compute_geometric_median_template(clips_s);
    ip=sum(sum(noise_shape.*subtemplate,1),2);
    subtemplate_resid=subtemplate-ip*noise_shape;
    subtemplate_resid_norm=sqrt(sum(subtemplate_resid(:).^2));
    subcluster.inds=inds_s;
    subcluster.detectability_score=subtemplate_resid_norm;
    subcluster.peak=subtemplate(channel,Tmid);
    subclusters{end+1}=subcluster;
end;
end


function noise_shape=estimate_noise_shape(X,T,ch)
% Computes the expected shape of the template in a noise cluster
% which may be considered as a row in the noise covariance matrix
% X is the MxN array of raw or pre-processed data
% T is the clip size
% ch is the channel where detection takes place
[M,N]=size(X);
clip_size=T;
Tmid=floor((T+1)/2);
num_rand_times=10000;
rand_times=clip_size+randsample(N-2*clip_size,num_rand_times)'; %random times
rand_clips=ms_extract_clips(X,rand_times,clip_size);
[~,~,L]=size(rand_clips);
peaks=rand_clips(ch,Tmid,:);

inds=find(abs(peaks)<=2); %focus on noise clips (where amplitude is low)
noise_shape=sum(rand_clips.*repmat(peaks,M,T,1),3);
noise_shape_norm=sqrt(sum(sum(noise_shape.^2,1),2));
noise_shape=noise_shape/noise_shape_norm;
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
