function mscmd_mda2txt(input_path,output_path,opts)
%MSCMD_COPY - Convert a .mda file to .txt file
%
% Syntax:  mscmd_mda2txt(input_path,output_path,opts)
%
% Inputs:
%    input_path - the path of the source .mda file
%    output_path - the path of the destination .txt file
%    opts -- script below
%
% Other m-files required: mscmd_exe
%
% See also: mscmd_copy.m

% Author: Jeremy Magland
% Mar 2016; Last revision: 4-Mar-2016

if (nargin<1) test_mscmd_mda2txt; return; end;
if (nargin<3) opts=struct; end;
if (~isfield(opts,'transpose')) opts.transpose=1; end;
if (~isfield(opts,'max_rows')) opts.max_rows=1e8; end;
if (~isfield(opts,'max_cols')) opts.max_cols=100; end;
if (~isfield(opts,'delim')) opts.delim='tab'; end; %can be 'comma'

cmd=sprintf('%s mda2txt --input=%s --output=%s --transpose=%d --max_rows=%ld --max_cols=%ld --delim=%s ',mscmd_exe,input_path,output_path,opts.transpose,opts.max_rows,opts.max_cols,opts.delim);

fprintf('\n*** MDA2TXT ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end

function test_mscmd_mda2txt
X=rand(4,50);
writemda(X,'tmp.mda');
opts.delim='comma';
mscmd_mda2txt('tmp.mda','tmp.txt',opts);
str=fileread('tmp.txt');
str
end

