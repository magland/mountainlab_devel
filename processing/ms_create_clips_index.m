function [clips_out,clips_index]=create_clips_index(clips,labels)

[M,T,NC]=size(clips);
clips_out=zeros(M,T,NC);
K=max(labels);
clips_index=zeros(1,K+1);
ind0=0;
for k=1:K
    inds_k=find(labels==k);
    clips_out(:,:,(ind0+1):(ind0+length(inds_k)))=clips(:,:,inds_k);
    clips_index(k)=ind0;
    ind0=ind0+length(inds_k);
end;
clips_index(end)=NC;

end
