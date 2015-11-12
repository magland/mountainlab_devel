function X=ms_whiten(X,opts)

M=size(X,1);
N=size(X,2);
num_clips=500;
clip_size=1;
npca=opts.num_whitening_components;
interval=ceil(N/num_clips);
times=1+clip_size:interval:N-clip_size;
clips=ms_extract_clips(X,times,clip_size);
[FF,info]=ms_event_features(clips,npca);

for cc=1:npca
    comp=squeeze(info.subspace(:,:,cc));
    ip=comp'*X; %(1xM) x (MxN) = 1xN
    X=X-comp*ip;
end;

end
