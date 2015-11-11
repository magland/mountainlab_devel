function view_cross_correlograms_example1

k=input('Enter a spike type number (or zero for diagonal): ');

output_path=example1_output_path;
path0=[output_path,'/cross-correlograms.mda'];
CC=mda_to_cross_correlograms(readmda(path0));
view_cross_correlograms(CC,k);

end

function CC=mda_to_cross_correlograms(X)

K=max(X(1,:));
CC=cell(K,K);

for k1=1:K
inds1=find(X(1,:)==k1);
for k2=1:K
    inds2=find(X(2,inds1)==k2);
    CC{k1,k2}=X(3,inds2);
end;
end;

end