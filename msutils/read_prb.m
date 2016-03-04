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

if nargin<1, test_read_dict; return; end;

str=fileread(fname);
X=parse_dict(str);

if nargout>=2
json=create_json(X);
end

end

function X=parse_dict(str)
tokens=tokenize_dict(str);
X=parse_dict_tokens(tokens);
end

function X=parse_dict_tokens(tokens)
if length(tokens)==0, X=[]; return; end;
if (length(tokens)==1),
    X=char(tokens{1});
    [~,status]=str2num(X); if status X=str2num(X); end;
    return; 
end;

if (strcmp(tokens{2},'='))
    subtokens=tokens(3:end);
    X=parse_dict_tokens(subtokens);
    return;
end;
if (strcmp(tokens{1},'['))||(strcmp(tokens{1},'('))
    if (~strcmp(tokens{end},closing_of(tokens{1})))
        error('Error parsing dict. Expected to find "%s".',closing_of(tokens{1}));
    end;
    append_entries=parse_delimiter_separated_array(tokens(2:end-1),'+');
    if (length(append_entries)==1)
        list=parse_delimiter_separated_array(append_entries{1},',');
        X1={};
        is_number_array=1;
        for kk=1:length(list)
            X0=parse_dict_tokens(list{kk});
            X1{end+1}=X0;
            if (~isnumeric(X0))||(length(X0)~=1)
                is_number_array=0;
            end;
        end;
        if is_number_array
            X=[];
            for j=1:length(X1)
                X(end+1)=X1{j};
            end;
        else
            X=X1;
        end;
    elseif (length(append_entries)>1)
        for kk=1:length(append_entries)
            X0=parse_dict_tokens(append_entries{kk});
            if (kk==1)
                if (isnumeric(X0)) X=[];
                else X={}; end;
            end;
            for j=1:length(X0)
                if (isnumeric(X0))
                    X(end+1)=X0(j);
                else
                    X{end+1}=X0{j};
                end;
            end;
        end;
    else
        error('Error parsing dict');
    end;
elseif (strcmp(tokens{1},'{'))
    if (~strcmp(tokens{end},closing_of(tokens{1})))
        error('Error parsing dict. Expected to find "%s".',closing_of(tokens{1}));
    end;
    entries=parse_delimiter_separated_array(tokens(2:end-1),',');
    X=struct;
    for j=1:length(entries)
        entry=entries{j};
        if (length(entry)==0) error('length of entry is zero'); end;
        if (strcmp(entry{1},'''')) %single quote, believe it or not
            if (length(entry)<5) error('length of entry should be >=5'); end;
            if (~strcmp(entry{3},'''')) error('unexpected dict parse problem'); end;
            if (~strcmp(entry{4},':')) error('unexpected dict parse problem'); end;
            subtokens=entry(5:end);
            X0=parse_dict_tokens(subtokens);
            name0=matlab.lang.makeValidName(char(entry{2}));
            X.(name0)=X0;
        else
            if (length(entry)<3) error('length of entry should be >=3'); end;
            if (~strcmp(entry{2},':')) error('unexpected dict parse problem'); end;
            subtokens=entry(3:end);
            X0=parse_dict_tokens(subtokens);
            name0=matlab.lang.makeValidName(char(entry{1}));
            X.(name0)=X0;
        end;
    end;
elseif (strcmp(tokens{1},'range'))
    X=[];
    if (length(tokens)<4) error('unexpected dict parse problem'); end;
    if (~strcmp(tokens{2},'(')) error('unexpected dict parse problem'); end;
    if (~strcmp(tokens{end},')'))
        error('unexpected dict parse problem');
    end;
    list=parse_delimiter_separated_array(tokens(3:end-1),',');
    if (length(list)==1)
        list1=list{1}; if (length(list1)~=1) error('unexpected dist parse problem'); end;
        num1=0; num2=str2num(char(list1{1}))-1;
    elseif (length(list)==2)
        list1=list{1}; if (length(list1)~=1) error('unexpected dist parse problem'); end;
        list2=list{2}; if (length(list2)~=1) error('unexpected dist parse problem'); end;
        num1=str2num(char(list1{1})); num2=str2num(char(list2{1}))-1;
    else
        error('unexpected dict parse problem');
    end
    for ii=num1:num2
        X(end+1)=ii;
    end;
else
    disp(tokens{1});
    disp(char(tokens{1}));
    error('Unexpected token.');
end;

end

function ret=parse_delimiter_separated_array(tokens,delim)
ret={};
level=0;
current_stuff={};
for j=1:length(tokens)
    if (is_begin_paren(tokens{j}))
        current_stuff{end+1}=tokens{j};
        level=level+1;
    elseif (is_end_paren(tokens{j}))
        current_stuff{end+1}=tokens{j};
        level=level-1;
        if (level<0)
            error('Unexpected dict parse problem');
        end;
    elseif (strcmp(tokens{j},delim))&&(level==0)
        if (length(current_stuff)>0)
            ret{end+1}=current_stuff;
            current_stuff={};
        else
            error('Unexpected dict parse problem.');
        end;
    else
        current_stuff{end+1}=tokens{j};
    end;
end;
if (level==0)
    if (length(current_stuff)>0)
        ret{end+1}=current_stuff;
        current_stuff={};
    end;
else
    error('Unexpected dict parse problem.');
end;
end

function ret=is_begin_paren(ch)
tmp='([{';
ret=~isempty(strfind(tmp,ch));
end

function ret=is_end_paren(ch)
tmp=')]}';
ret=~isempty(strfind(tmp,ch));
end

function ret=closing_of(ch)
if (strcmp(ch,'[')) ret=']'; return; end;
if (strcmp(ch,'(')) ret=')'; return; end;
if (strcmp(ch,'{')) ret='}'; return; end;
error('Unexpected character in closing_of: "%s"',ch);
end

function tokens=tokenize_dict(str)
tmp='(\-\d+\.\d+|\-\d+|\d+\.\d+|\w+|\+|{|}|\(|\)|\[|\]|\:|\=|,|\-|'')';
tokens=regexp(str,tmp,'tokens');
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

function test_read_dict
mfile_path=fileparts(mfilename('fullpath'));
[X,json]=read_prb([mfile_path,'/sample.prb']);
json
X
X_from_json=JSON.parse(json)
end