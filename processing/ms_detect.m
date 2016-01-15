function [times_pos,times_neg]=ms_detect(X,opts)

if nargin==0 test_ms_detect; return; end;

[times_pos]=ms_detect_2(X,opts);
times_neg=[];
return;

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

function [times]=ms_detect_2(X,opts)

detect_interval=opts.detect_interval;
detect_threshold=opts.detect_threshold;

absX=abs(X);
absX=max(absX,[],1); %max over channels (but best to put in one channel at a time)

N=length(X);
use_it=zeros(1,N);
best_ind=1;
best_abs_val=absX(1);
candidates=find((absX>=detect_threshold)&((1:N)>opts.clip_size)&((1:N)<=N-opts.clip_size));
for tt=candidates
    if (best_ind<tt-detect_interval)
        [~,best_ind]=max(absX(tt-detect_interval:tt-1));
        best_ind=best_ind+tt-detect_interval-1;
        best_abs_val=absX(best_ind);
    end;
    if (absX(tt)>best_abs_val)
        use_it(tt)=1;
        use_it(best_ind)=0;
        best_ind=tt;
        best_abs_val=absX(tt);
    end;
end;

times=find(use_it==1);

end

function test_ms_detect
opts.detect_threshold=10;
opts.detect_interval=5;
opts.clip_size=50;
X=[rand(1,150)*10,rand(1,150)*20];
times=ms_detect(X,opts);
figure; plot(1:length(X),X,'k'); hold on;
for j=1:length(times)
    plot([times(j),times(j)],ylim,'r');
end;
end
