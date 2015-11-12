function [CC,CCmda]=ms_cross_correlograms(times,labels,max_dt)

[times,inds]=sort(times);
labels=labels(inds);

K=max(labels);
CC=cell(K,K);

i1=1;
for i2=1:length(times)
    while (times(i1)<times(i2)-max_dt) i1=i1+1; end;
    k2=labels(i2);
    t2=times(i2);
    for jj=i1:i2-1
        k1=labels(jj);
        t1=times(jj);
        CC{k1,k2}=[CC{k1,k2},t2-t1];
        CC{k2,k1}=[CC{k2,k1},t1-t2];
    end;
end;

CCmda=cross_correlograms_to_mda(CC);

end

function ret=cross_correlograms_to_mda(CC)
K=size(CC,1);

ct=0;
for k1=1:K
for k2=1:K
ct=ct+length(CC{k1,k2});
end;
end;
ret=zeros(3,ct);

ct=0;
for k1=1:K
for k2=1:K
ret(1,ct+1:ct+length(CC{k1,k2}))=k1;
ret(2,ct+1:ct+length(CC{k1,k2}))=k2;
ret(3,ct+1:ct+length(CC{k1,k2}))=CC{k1,k2};
ct=ct+length(CC{k1,k2});
end;
end;

end

