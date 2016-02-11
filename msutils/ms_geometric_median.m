function [mm,changes]=ms_geometric_median(X,num_iterations)
if nargin<2, num_iterations=10; end;

[M,N]=size(X);
weights=ones(1,N);
changes=[];
for it=1:num_iterations
    weights=weights/sum(weights);
    mm=X*weights';
    if (it>1)
        changes=[changes,sqrt(sum((mm-mm_old).^2))];
    end;
    mm_old=mm;
    diffs=X-repmat(mm,1,N);
    weights=sqrt(sum(diffs.^2,1));
    inds=find(weights~=0);
    weights(inds)=1./weights(inds);
end;

end

