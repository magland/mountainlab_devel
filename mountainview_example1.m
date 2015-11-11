function mountainview_example1

output_path=[fileparts(mfilename('fullpath')),'/../ms11d45A/output'];
%mountainview(output_path); return;

times=readmda([output_path,'/times.mda']);
labels=readmda([output_path,'/labels.mda']);
CC=compute_cross_correlograms(times,labels,1500);

figure;
for j=1:7*7
    subplot(7,7,j);
    if (j<=size(CC,2))
        hist(CC{43,j},250);
        title(sprintf('%d',j));
    end;
end;

view_cross_correlograms(CC,43);

CCmda=readmda([output_path,'/cross-correlograms.mda']);
CC2=mda_to_cross_correlograms(CCmda);
view_cross_correlograms(CC2,43);

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