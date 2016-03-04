function [X,json]=read_prm(fname)
%READ_PRM - read a python dict file, returning a matlab structure,
%cell array, or numerical array. This is probably VERY buggy, but should
%read .prm files for spike sorting
%
% Syntax:  [X] = read_prm(fname)
%
% Inputs:
%    fname - the file name for the python dict, e.g., *.prm
%
% Outputs:
%    X - the matlab struct, cell array, or numerical arrays
%    json - JSON text of X
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
% See also: read_prb.m

% Author: Jeremy Magland
% Mar 2015; Last revision: 4-Mar-2016

if nargin<1, test_read_prm; return; end;

str=fileread(fname);
X=parse_prm(str);

if nargout>=2
json=create_json(X);
end
end

function X=parse_prm(str)
lines=strsplit(str,sprintf('\n'));
str2='{';
for j=1:length(lines)
    line=lines{j};
    ind=strfind(line,'#');
    if (length(ind)>=1)
        line=line(1:ind(1)-1);
    end;
    line=strtrim(line);
    if (~isempty(strfind(line,'*')))||(~isempty(strfind(line,'+')))
        %do not allow in-line arithmetic!
        line='';
    end;
    if (length(line)>0)
        
        str2=[str2,line];
        if (~strcmp(line(end),','))&&(~strcmp(line(end),'('))
            str2=[str2,','];
        end;
        str2=[str2,sprintf('\n')];
    end;
end;
str2=[str2,'}'];
str2=strrep(str2,'=',':');
str2=strrep(str2,'=',':');
str2=strrep(str2,'dict(','{');
str2=strrep(str2,')','}');
X=parse_prb(str2);
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

function test_read_prm
mfile_path=fileparts(mfilename('fullpath'));
[X,json]=read_prm([mfile_path,'/sample.prm']);
json
X
X_from_json=JSON.parse(json)
end