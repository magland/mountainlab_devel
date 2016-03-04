function mda2txtfile(mda_fname,txt_fname,opts)

if nargin<1 test_mda2txtfile; return; end;

if nargin<3, opts=struct; end;
if (~isfield(opts,'delimiter')) opts.delimiter=sprintf('\t'); end;
if (~isfield(opts,'transpose')) opts.transpose=1; end;
if (~isfield(opts,'max_rows')) opts.max_rows=1e9; end;
if (~isfield(opts,'max_cols')) opts.max_cols=100; end;

X=readmda(mda_fname);
if ndims(X)>2
error('Number of dimensions is larger than 2');
end;
if opts.transpose X=X'; end;

if (size(X,1)>opts.max_rows)
warning('Unable to create mda text file, too many rows.');
return;
end;

if (size(X,2)>opts.max_cols)
warning('Unable to create mda text file, too many columns.');
return;
end;

dlmwrite(txt_fname,X,opts.delimiter);
end

function test_mda2txtfile
X=rand(4,50);
writemda(X,'tmp.mda');
mda2txtfile('tmp.mda','tmp.txt');
str=fileread('tmp.txt');
str
end