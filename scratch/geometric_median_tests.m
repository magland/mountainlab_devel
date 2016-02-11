function geometric_median_tests

close all;

N1=500;
N2=100;
N3=30;
X1=randn(2,N1);
X2=randn(2,N2) + repmat([30,15]',1,N2);
X3=randn(2,N3) + repmat([-20,12]',1,N3);
X=cat(2,X1,X2,X3);
[mm,changes]=geometric_median(X,10);
figure; semilogy(1:length(changes),changes);

cc=mean(X,2);

figure;
ms_view_clusters(X); hold on;
plot([mm(1),mm(1)],[mm(2),mm(2)],'r.','MarkerSize',20); %geometric median
plot([cc(1),cc(1)],[cc(2),cc(2)],'g.','MarkerSize',20); %centroid

legend('data','geometric median','centroid');

end

function [mm,changes]=geometric_median(X,num_iterations)

if nargin<2, num_iterations=10; end;

[M,N]=size(X);
weights=rand(1,N);
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