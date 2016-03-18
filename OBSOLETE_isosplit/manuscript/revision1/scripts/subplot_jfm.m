function subplot_jfm(nr,nc,ii,opts)

if (nargin<4) opts=struct; end;
if (~isfield(opts,'xmargin')) opts.xmargin=0.01; end;
if (~isfield(opts,'ymargin')) opts.ymargin=0.05; end;
if (~isfield(opts,'xspacing')) opts.xspacing=0.02; end;
if (~isfield(opts,'yspacing')) opts.yspacing=0.05; end;
if (~isfield(opts,'yoffset')) opts.yoffset=0; end;

set(gcf,'Color','w');

c0=mod((ii-1),nc)+1;
r0=floor((ii-1)/nc)+1;
r0=nr+1-r0;

xspacing=opts.xspacing; yspacing=opts.yspacing;
yoffset=opts.yoffset;
xmargin=opts.xmargin; ymargin=opts.ymargin;
W0=(1-(xspacing*(nc-1)+xmargin*2))/nc;
H0=(1-(yspacing*(nr-1)+ymargin*2))/nr;

subplot('Position',[xmargin+(c0-1)*(W0+xspacing),ymargin+(r0-1)*(H0+yspacing)+yoffset,W0,H0]);

end
