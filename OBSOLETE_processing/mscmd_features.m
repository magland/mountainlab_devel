function mscmd_features(input_clips_path,detect_path,adjacency_path,output_path,opts)
%MSCMD_FEATURES - Extract PCA features from an array of events
%
%This is a wrapper to the command-line mountainsort command. It has
%similar functionality to ms_event_features.
%
% Syntax:  mscmd_features(input_clips_path,detect_path,adjacency_path,output_path,opts)
%
% Inputs:
%    input_clips_path - path to MxTxNC array of clips (see mscmd_extract_clips)
%    detect_path - 2xNC array output from mscmd_detect containing the
%                  event channel numbers in the first row
%    adjacency_path - path to the MxM adjacency matrix of 0's and 1's
%                     defining which channels are adjacent to one another.
%                     The PCA analysis will only be performed on a channel
%                     and its adjacent neighbors.
%    output_path - path to the (num_features x NC) array of features
%    opts.num_features - the number of features to extract
%    opts.clip_size - why is this parameter needed?
%
% Other m-files required: mscmd_exe
%
% See also: ms_event_features, mscmd_extract_clips, mscmd_detect

% Author: Jeremy Magland and Alex Barnett
% Oct 2015; Last revision: 13-Feb-2016

if (nargin<4) opts=struct; end;

cmd=sprintf('%s features --input=%s --detect=%s --adjacency=%s --output=%s ',mscmd_exe,input_clips_path,detect_path,adjacency_path,output_path);
cmd=[cmd,sprintf('--num_features=%d --clip_size=%d ',opts.num_features,opts.clip_size)];

fprintf('\n*** FEATURES ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end