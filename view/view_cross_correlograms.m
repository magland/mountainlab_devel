function view_cross_correlograms(CC,k)

if (isstr(CC))
    CC=readmda(CC);
end;
if (~iscell(CC))
    CC=mda_to_cross_correlograms(CC);
end;

K=size(CC,1);
aa=ceil(sqrt(K));

figure('name','Cross-correlograms');
set(gcf,'position',[100,100,1000,1000]);
for j=1:aa*aa
    k0=k;
    if (k0==0) k0=j; end;
    if (j<=size(CC,2))
        subplot(aa,aa,j);
        hist(CC{k0,j},250);
        count=length(CC{k0,j});
        title(sprintf('%d / %d (%d)',k0,j,count),'fontsize',8);
        set(gca,'xticklabel',[]);
        set(gca,'yticklabel',[]);       
    end;
end;

end

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
