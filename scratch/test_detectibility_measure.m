function test_detectibility_measure

close all;

clip_size=100;

path0='franklab/test_hippocampal_02_10_2016/output_tetrode1';
pre2=readmda([path0,'/pre2.mda']);
firings=readmda([path0,'/firings.mda']);
times=firings(2,:);
labels=firings(3,:);
peaks=firings(4,:);
%inds=find((abs(peaks)<6)&(abs(peaks)>5));
inds=find(abs(peaks)>0);
times=times(inds); labels=labels(inds); peaks=peaks(inds);
clips=ms_extract_clips(pre2,times,clip_size);
templates=ms_templates(clips,labels);
figure; ms_view_templates(templates);
K=size(templates,3);

detectibility_scores=zeros(1,K);
template_norms=zeros(1,K);

for k=2
template0=templates(:,:,k);
%figure; ms_view_templates(template0);

num_rclips=5000;
rclips=extract_random_clips(pre2,clip_size,num_rclips);
inner_products=squeeze(sum(sum(repmat(template0,1,1,num_rclips).*rclips,1),2));
sigma=sqrt(var(inner_products));
normsqr_template0=sum(template0(:).^2);
figure; hist(inner_products,1000);
hold on;
plot([normsqr_template0,normsqr_template0],ylim,'r','linewidth',15);
detectibility_scores(k)=normsqr_template0/sigma
template_norms(k)=sqrt(normsqr_template0)
end;

figure; plot(1:K,template_norms,'b.',1:K,detectibility_scores,'r.','markersize',8);
figure; plot(template_norms,detectibility_scores,'k.','markersize',8);

end

function clips=extract_random_clips(X,clip_size,num_clips)
[M,N]=size(X);
times=randi([clip_size+1,N-clip_size-1],1,num_clips);
clips=ms_extract_clips(X,times,clip_size);
end
