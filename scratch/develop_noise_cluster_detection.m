function develop_noise_cluster_detection

close all;

template1=test_1([0,inf]);
figure; ms_view_templates(template1);

template2=test_1([500,3000]);
figure; ms_view_templates(template2);

template3=test_2([500,3000]);
figure; ms_view_templates(cat(3,template2,template3));

[templates,templates_resid,detectability_scores,peaks]=test_3;
figure; ms_view_templates(templates); 
title('Cluster Templates');
figure; ms_view_templates(templates_resid); 
title('Residual templates after subtracting estimated noise shapes');
figure; bar(1:length(detectability_scores),detectability_scores);
xlabel('Cluster'); ylabel('Detectability Score'); title('Detectability Scores');
figure; plot(abs(peaks),detectability_scores,'b.','MarkerSize',8);
xlabel('Abs. Peak Amplitude'); ylabel('Detectability Score');


end

function [template,X]=test_1(freq_range)

M=4;
N=1e5;
clip_size=100; T=clip_size; Tmid=floor((T+1)/2);
ch=1;
X=randn(M,N);

o_filter.samplefreq=20000;
o_filter.freq_min=freq_range(1);
o_filter.freq_max=freq_range(2);
X=ms_filter(X,o_filter);
X=ms_whiten(X);

o_detect.detect_threshold=3;
o_detect.detect_interval=50;
o_detect.clip_size=clip_size;
times=ms_detect(X(ch,:),o_detect);
inds_pos=find(X(ch,times)>0); times=times(inds_pos);
clips=ms_extract_clips(X,times,clip_size);
[~,~,L]=size(clips);
peaks=reshape(clips(ch,Tmid,:),1,L);

template=mean(clips,3);

end

function [template,X,noise_shape]=test_2(freq_range)

ch=1;

[template,X]=test_1(freq_range);
[M,T]=size(template); [M,N]=size(X);
noise_shape=estimate_noise_shape(X,T,ch);
ip=sum(sum(noise_shape.*template,1),2);
template=template-ip*noise_shape;

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

function [templates,templates_resid,detectability_scores,peaks]=test_3
mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../franklab/sort_003_multichannel/output_tetrode1'];

clip_size=100;

pre2=readmda([path0,'/pre2.mda']);
[M,N]=size(pre2);
T=clip_size; Tmid=floor((T+1)/2);
firings=readmda([path0,'/firings.mda']);
channels=firings(1,:); times=firings(2,:); labels=firings(3,:); peaks=firings(4,:);

clips=ms_extract_clips(pre2,times,clip_size);
templates=ms_templates(clips,labels);
K=size(templates,3);
figure; ms_view_templates(templates);

templates_resid=zeros(size(templates));
peaks=zeros(1,K);
for k=1:K
    disp(k);
    inds_k=find(labels==k);
    ch=channels(inds_k(1));
    noise_shape=estimate_noise_shape(pre2,clip_size,ch);
    template=templates(:,:,k);
    ip=sum(sum(noise_shape.*template,1),2);
    template_resid=template-noise_shape*ip;
    templates_resid(:,:,k)=template_resid;
    peaks(k)=template(ch,Tmid,:);
end;

resids=reshape(sqrt(sum(sum(templates_resid.^2,1),2)),1,K);
detectability_scores=resids;
coeffs=polyfit(abs(peaks),detectability_scores,1);
detectability_scores=detectability_scores/coeffs(1);

end
