function CC=mda_to_cross_correlograms(X)

K=max(X(1,:));
CC=cell(K,K);

for k1=1:K
inds1=find(X(1,:)==k1);
for k2=1:K
    inds2=find(X(2,inds1)==k2);
    CC{k1,k2}=X(3,inds1(inds2));
end;
end;

end 
