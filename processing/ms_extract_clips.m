function clips=extract_clips(Y,times,clip_size)
% Extract clips centered at times

[M,N]=size(Y);
T=clip_size;
C=length(times);

clips=zeros(M,T,C);
tt1=-ceil((clip_size)/2);
tt2=tt1+clip_size-1;
inds=find((times+tt1>=1)&(times+tt2<=N));
%if (min(times+tt1)<1) error('Invalid time in extract_clips'); end;
%if (max(times+tt2)>N) error('Invalid time in extract_clips'); end;
%for j=1:C
for j=inds
	clips(:,:,j)=Y(:,times(j)+tt1:times(j)+tt2);
end;

end
