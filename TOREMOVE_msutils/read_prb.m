function [X,json]=read_prb(fname)
%READ_PRB - read a python dict file, returning a matlab structure,
%cell array, or numerical array. This is probably VERY buggy, but should
%read .prb files for spike sorting
%
% Syntax:  [X] = read_prb(fname)
%
% Inputs:
%    fname - the file name for the python dict, e.g., *.prb
%
% Outputs:
%    X - the matlab struct, cell array, or numerical arrays
%
% The fields of matlab structs must not be numbers, so these are converted
% to valid field names by prepending the 'x' characters, e.g., X.x0 rather
% than X.0. More specifically we use matlab.lang.makeValidName.
%
% For convenience, cell arrays are converted to numerical arrays whenever
% possible.
%
% Other m-files required: none
%
% See also: 

% Author: Jeremy Magland
% Mar 2015; Last revision: 4-Mar-2016

if nargin<1, test_read_prb; return; end;

str=fileread(fname);
X=parse_prb(str);

if nargout>=2
json=create_json(X);
end

end


function json=create_json(X)
json='';
if (isstruct(X))
    json=[json,'{'];
    FF=fieldnames(X);
    for j=1:length(FF)
        tmp=create_json(X.(FF{j}));
        json=[json,sprintf('"%s":%s',FF{j},tmp)];
        if (j+1<=length(FF)) json=[json,',']; end;
    end
    json=[json,'}'];
elseif (iscell(X))
    json=[json,'['];
    for j=1:length(X)
        tmp=create_json(X{j});
        json=[json,tmp];
        if (j+1<=length(X)) json=[json,',']; end;
    end
    json=[json,']'];
elseif (isnumeric(X))&&(length(X)==1)
    json=[json,num2str(X)];
elseif (isnumeric(X))
    json=[json,'['];
    for j=1:length(X)
        json=[json,num2str(X(j))];
        if (j+1<=length(X)) json=[json,',']; end;
    end
    json=[json,']'];
elseif (isstr(X))
    json=[json,'"',X,'"'];
else
    error('Unexpected type.');
end
end

function test_read_prb
mfile_path=fileparts(mfilename('fullpath'));
[X,json]=read_prb([mfile_path,'/sample.prb']);
json
X
X_from_json=JSON.parse(json)
end