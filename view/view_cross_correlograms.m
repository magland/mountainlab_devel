function view_cross_correlograms(CC,k)
%VIEW_CROSS_CORRELOGRAMS - View cross-correlograms
%
% Syntax:  view_cross_correlograms(CC,k)
%
% Inputs:
%    CC - either CC or CCmda output from ms_cross_correlograms, or CCmda
%    output from mscmd_cross_correlograms
%    k - either zero (for auto-correlograms) or >0 for cross-correlograms
%    with label k
%
% Other m-files required: none
%
% See also: ms_cross_correlograms, mscmd_cross_correlograms,
% ms_mountainview

% Author: Jeremy Magland
% Jan 2015; Last revision: 15-Feb-2106

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

