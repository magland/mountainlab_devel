function [times_pos,times_neg]=ms_detect(X,opts)

detect_interval=opts.detect_interval;
detect_threshold=opts.detect_threshold;

absX=abs(X);
absX=max(absX,[],1); %max over channels (but best to put in one channel at a time)
times=find(absX>detect_threshold);
signed_vals=X(times);
pos_inds=find(signed_vals>=0);
neg_inds=find(signed_vals<0);

%Now make sure we only use those that are global maxima over the radius
%Nt/2
use_it=ones(1,length(times));
for j=1:length(times)
    if ((times(j)-ceil(opts.clip_size/2))<1) use_it(j)=0; end;
    if ((times(j)+ceil(opts.clip_size/2))>size(X,2)) use_it(j)=0; end;
    k=j-1;
    while (k>=1)&&(times(k)>=times(j)-detect_interval)
        if (absX(times(k))>absX(times(j)))
            use_it(j)=0;
        else
            use_it(k)=0;
        end;
        k=k-1;
    end
end;

times_pos=times(pos_inds);
times_neg=times(neg_inds);
times_pos=times_pos(find(use_it(pos_inds)));
times_neg=times_neg(find(use_it(neg_inds)));

end
