function test_detect_2_26_2016

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../franklab/sort_003_multichannel/output_tetrode1'];

clip_size=200;

pre2=readmda([path0,'/pre2.mda']);
[M,N]=size(pre2);
T=clip_size; Tmid=floor((T+1)/2);
firings=readmda([path0,'/firings.mda']);
channels=firings(1,:); times=firings(2,:); labels=firings(3,:); peaks=firings(4,:);

clips=ms_extract_clips(pre2,times,clip_size);
templates=ms_templates(clips,labels);
figure; ms_view_templates(templates);

num_rand_times=50000;
rand_times=clip_size+randsample(N-2*clip_size,num_rand_times,true)';
rand_clips=ms_extract_clips(pre2,rand_times,clip_size);

k=4;
template=templates(:,:,k);
inds_k=find(labels==k);
ch=channels(inds_k(1));
template_peak=template(ch,Tmid);
rand_peaks=reshape(rand_clips(ch,Tmid,:),1,num_rand_times);
inds00=find((0<rand_peaks*sign(template_peak))&(rand_peaks*sign(template_peak)<2));
shape=mean(rand_clips(:,:,inds00),3);
norm_shape=sqrt(sum(sum(shape.^2,1),2));
shape=shape/norm_shape;
ip=sum(sum(shape.*template,1),2);
figure; ms_view_templates(shape);
figure; ms_view_templates(cat(3,template,template-shape*ip));

return;

num_rand_times=50000;
rand_times=clip_size+randsample(N-2*clip_size,num_rand_times,true)';
rand_clips=ms_extract_clips(pre2,rand_times,clip_size);
rand_peaks=reshape(rand_clips(ch,Tmid,:),1,num_rand_times);
inds=find((0<rand_peaks)&(rand_peaks<2));
shape=mean(rand_clips(:,:,inds),3);
figure; ms_view_templates(shape);
Vshape=reshape(shape,M*T,1);
Vshape_norm=sqrt(Vshape'*Vshape);
Vshape=Vshape/Vshape_norm;

times=ms_detect(pre2(ch,:),o_detect);
L=length(times);
clips=ms_extract_clips(pre2,times,clip_size);
[M,T,L]=size(clips); Tmid=floor((T+1)/2);
Vclips=reshape(clips,M*T,L);
ips=Vshape'*Vclips;
clipsB=clips-repmat(reshape(ips,1,1,L),M,T,1).*repmat(shape/Vshape_norm,1,1,L);

FF=ms_event_features(clips,3);
labels=isosplit2(FF);
figure; ms_view_clusters(FF,labels);
figure; ms_view_templates_from_clips(clips,labels);

FFB=ms_event_features(clipsB,3);
labelsB=isosplit2(FFB);
figure; ms_view_clusters(FFB,labelsB);
figure; ms_view_templates_from_clips(clips,labelsB);


return;



times=ms_detect(pre2(ch,:),o_detect);
clips=ms_extract_clips(pre2,times,clip_size);
[M,T,L]=size(clips); Tmid=floor((T+1)/2);
peaks=reshape(clips(ch,Tmid,:),1,L);

indsA=find((2<=peaks)&(peaks<3));
indsB=find((4<=peaks));
clipsA=clips(:,:,indsA);
clipsB=clips(:,:,indsB);
[~,~,LB]=size(clipsB);
noise_peak_shape=mean(clipsA,3);
V0=reshape(noise_peak_shape,M*T,1);
norm0=sqrt(V0'*V0);
V0=V0/norm0;
VclipsB=reshape(clipsB,M*T,LB);
ips=V0'*VclipsB;
clipsB2=clipsB-repmat(reshape(ips,1,1,LB),M,T,1).*repmat(noise_peak_shape/norm0,1,1,LB);
tmp=mean(clipsB2(:,:,:),3);
figure; ms_view_templates(cat(3,noise_peak_shape,tmp));

FF1=ms_event_features(clipsB,3);
labels1=isosplit2(FF1);
figure; ms_view_clusters(FF1,labels1);

FF2=ms_event_features(clipsB2,3);
labels2=isosplit2(FF2);
figure; ms_view_clusters(FF2,labels2);

end
