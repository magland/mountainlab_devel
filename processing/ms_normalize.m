function X=ms_normalize(X)
%MS_NORMALIZE - Normalize the data so each channel has unit variance
%
%Consider using ms_whiten or mscmd_whiten
%
% Syntax:  [Y] = ms_normalize(X)
%
% Inputs:
%    X - MxN array of raw or pre-processed data
%
% Outputs:
%    Y - MxN array of normalized data
%
% Other m-files required: none
%
% See also: ms_whiten, mscmd_whiten

% Author: Jeremy Magland
% Jan 2016; Last revision: 13-Feb-2016
for j=1:size(X,1)
    stdev=sqrt(var(X(j,:)));
    X(j,:)=X(j,:)/stdev;
    X(j,:)=min(20,max(-20,X(j,:)));
end;
for j=1:size(X,1)
    stdev=sqrt(var(X(j,:)));
    X(j,:)=X(j,:)/stdev;
end;
end