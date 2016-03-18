function mscmd_normalize_channels(input_path,output_path)
%MSCMD_NORMALIZE_CHANNELS - Normalize the data so each channel has unit variance
%
%Consider using mscmd_whiten. This is a wrapper to the command-line
%mountainsort command with functionality similar to ms_normalize
%
% Syntax:  mscmd_normalize_channels(input_path,output_path)
%
% Inputs:
%    input_path - path to MxN array of raw or pre-processed data
%    output_path - path to MxN normalized array
%
% See also: mscmd_whiten, ms_normalize

% Author: Jeremy Magland
% Jan 2016; Last revision: 13-Feb-2016

cmd=sprintf('%s normalize_channels --input=%s --output=%s ',mscmd_exe,input_path,output_path);

fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end