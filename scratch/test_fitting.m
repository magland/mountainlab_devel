function test_fitting

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0='franklab/test_hippocampal_01_28_2016/tetrode1_output';

clips0=readmda([path0,'/clips0.mda']);
clips=readmda([path0,'/clips.mda']);
[M,T,NC]=size(clips);
clips_index=readmda([path0,'/clips_index.mda']);
templates=readmda([path0,'/templates.mda']);
size(templates)
%ms_view_templates(templates);
k=10;
inds=clips_index(k)+1:clips_index(k+1);
%template0=templates(:,:,k);
template0=clips(:,:,randi(size(clips,3)));
template0=template0/sqrt(sum(template0(:).^2));
ips=squeeze(sum(sum(repmat(template0,1,1,NC).*clips,1),2));
ips0=squeeze(sum(sum(repmat(template0,1,1,size(clips0,3)).*clips0,1),2));
norms=squeeze(sqrt(sum(sum(clips.^2,1),2)));
norms0=squeeze(sqrt(sum(sum(clips0.^2,1),2)));
figure; 
plot(norms0,ips0,'b.'); hold on;
%plot(norms(inds),ips(inds),'r.'); 
figure; ms_view_templates(template0);

%tmp=cat(1,norms0',ips0');
%labels0=isosplit(tmp,struct('K',100));
%ms_view_clusters(tmp,labels0);


return;




clips=clips(:,:,clips_index(4):clips_index(5)-1);
spikespy(clips);
%size(clips)
aa=clips([1,3],101,:);
%figure; hist(squeeze(aa(1,:,:)),1000);
%figure; hist(squeeze(aa(2,:,:)),1000);
figure; plot(aa(1,:),aa(2,:),'b.');
return;
FF=ms_event_features(aa,6);
ms_view_clusters(FF);
return;

if 0
    pre2=readmda([path0,'/pre2.mda']);
    clusters=readmda([path0,'/clusters.mda']);

    pre2=pre2(:,1:1e5);
    inds=find(clusters(2,:)<=size(pre2,2));
    clusters=clusters(:,inds);
    writemda(pre2,[mfile_path,'/test_fitting_pre2.mda']);
    writemda(clusters,[mfile_path,'/clusters.mda']);
end;

pre2=readmda([mfile_path,'/test_fitting_pre2.mda']);
clusters=readmda([mfile_path,'/clusters.mda']);
templates=readmda([path0,'/templates.mda']);

times=clusters(2,:);
labels=clusters(3,:);

figure;
ms_view_templates(templates);

%figure;
%hist(labels,1:max(labels));

k=4;
inds_k=find(labels==k);
aa=compute_sliding_ip(pre2,templates);
for m=1:size(aa,1)
    aa(m,:)=aa(m,:)/sqrt(var(aa(m,:)));
end;

%spikespy({aa,times(inds_k),labels(inds_k)});
%spikespy({cat(1,pre2,aa),times(inds_k),labels(inds_k)});

%spikespy({aa,times,labels});

times0=ms_detect(aa,struct('detect_threshold',4,'detect_interval',15,'clip_size',200));
clips0=ms_extract_clips(aa,times0,100);
FF=ms_event_features(clips0,6);
ms_view_clusters(FF);

end

function ret=compute_sliding_ip(X,kernel)
[M,N]=size(X);
[~,T,K]=size(kernel);

if K>1
    ret=zeros(K,N);
    for k=1:K
        ret(k,:)=compute_sliding_ip(X,kernel(:,:,k));
    end;
    return;
end;

tt1=-T/2;
tt2=tt1+T-1;

ret=zeros(M,N);
for shift=tt1:tt2
    if (shift>=0)
        ret(:,1:N-shift)=ret(:,1:N-shift)+repmat(kernel(:,shift-tt1+1),1,N-shift).*X(:,1+shift:N);
    else
        ret(:,1-shift:N)=ret(:,1-shift:N)+repmat(kernel(:,shift-tt1+1),1,N+shift).*X(:,1:N+shift);
    end;
end;

ret=sum(ret,1);

end

