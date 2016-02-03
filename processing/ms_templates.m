function templates=ms_templates(clips,labels)

if (nargin<3) opts=struct; end;

[M,T,NC]=size(clips);

K=max(labels);

templates=zeros(M,T,K);

for k=1:K
    ii0=find(labels==k);
    clips0=clips(:,:,ii0);
    templates(:,:,k)=mean(clips0,3);
end;

end