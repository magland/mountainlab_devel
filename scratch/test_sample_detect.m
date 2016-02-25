function test_sample_detect

N=800;
detect_interval=15;
detect_threshold=4;
X=randn(1,N)+create_triangular_spikes(N,10,5,30);
times=sample_detect(X,detect_threshold,detect_interval);
figure; plot(1:N,X,'k'); hold on;
for j=1:length(times)
    plot([times(j),times(j)],ylim,'r--');
end;
plot(xlim,[detect_threshold,detect_threshold],'b--');
xlabel('Timepoints');
ylabel('Signal');

end

function times=sample_detect(X,detect_threshold,detect_interval)
N=length(X);
use_it=zeros(1,N);
best_ind=1;
best_val=X(1);
candidates=find(X>=detect_threshold);
for tt=candidates
    if (best_ind<tt-detect_interval)
        [~,best_ind]=max(X(tt-detect_interval:tt-1));
        best_ind=best_ind+tt-detect_interval-1;
        best_val=X(best_ind);
    end;
    if (X(tt)>=best_val)
        use_it(tt)=1;
        use_it(best_ind)=0;
        best_ind=tt;
        best_val=X(tt);
    end;
end;

times=find(use_it==1);
end

function X=create_triangular_spikes(N,num,amp,len)
X=zeros(1,N);
spike_shape=amp*(1-abs(linspace(-1,1,len)));
for n=1:num
    loc=randi([1,N-len+1]);
    X(loc:loc+len-1)=spike_shape;
end;
end