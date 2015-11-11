function X=compute_cross_correlograms(times,labels,max_dt)

[times,inds]=sort(times);
labels=labels(inds);

K=max(labels);
X=cell(K,K);

i1=1;
for i2=1:length(times)
    while (times(i1)<times(i2)-max_dt) i1=i1+1; end;
    k2=labels(i2);
    t2=times(i2);
    for jj=i1:i2-1
        k1=labels(jj);
        t1=times(jj);
        X{k1,k2}=[X{k1,k2},t2-t1];
        X{k2,k1}=[X{k2,k1},t1-t2];
    end;
end;

end
