function CC=mda_to_cross_correlograms(CCmda)
%MDA_TO_CROSS_CORRELOGRAMS - utility to convert from CCmda format to CC
%(see ms_cross_correlograms and mscmd_cross_correlograms)
%
% Syntax:  [CC] = mda_to_cross_correlograms(X)
%
% Inputs:
%    CCmda - output from ms_cross_correlograms, or output from 
%        mscmd_cross_correlograms
%
% Outputs:
%    CC - cell array containing cross-correlogram data
%
% Other m-files required: none
%
% See also: ms_cross_correlograms, mscmd_cross_correlograms,
% ms_mountainview

% Author: Jeremy Magland
% Jan 2015; Last revision: 15-Feb-2106

K=max(CCmda(1,:));
CC=cell(K,K);

for k1=1:K
inds1=find(CCmda(1,:)==k1);
for k2=1:K
    inds2=find(CCmda(2,inds1)==k2);
    CC{k1,k2}=CCmda(3,inds1(inds2));
end;
end;

end 
