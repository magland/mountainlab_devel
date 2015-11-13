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

