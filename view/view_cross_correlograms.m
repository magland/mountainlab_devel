function view_cross_correlograms(CC,k)

K=size(CC,1);
aa=ceil(sqrt(K));

figure('name','Cross-correlograms');
for j=1:aa*aa
    k0=k;
    if (k0==0) k0=j; end;
    subplot(aa,aa,j);
    if (j<=size(CC,2))
        hist(CC{k0,j},250);
        %title(sprintf('%d / %d',k0,j));
    end;
end;

end
