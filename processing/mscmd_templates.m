function mscmd_templates(input_clips_path,clusters_path,output_path,opts)
%MSCMD_TEMPLATES - Compute templates (mean waveforms) corresponding to each
%spike type.
%
%This is a wrapper to the command-line mountainsort call, with
%functionality similar to ms_templates.
%
% Syntax:  mscmd_templates(input_clips_path,clusters_path,output_path,opts)
%
% Inputs:
%    input_clips_path - path to MxTxNC array of clips
%    clusters_path - path to array containing the times/labels.
%                    See docs for format of this array
%    output_path - path to MxTxK array of computed templates, K is the
%                  number of spike types (max label)
%
% Other m-files required: mscmd_exe
%
% See also: ms_templates, mscmd_extract_clips

% Author: Jeremy Magland
% Jan 2015; Last revision: 13-Feb-2106

if (nargin<4) opts=struct; end;

cmd=sprintf('%s templates --input=%s --clusters=%s --output=%s --clip_size=%d ',mscmd_exe,input_clips_path,clusters_path,output_path,opts.clip_size);

fprintf('\n*** TEMPLATES ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end