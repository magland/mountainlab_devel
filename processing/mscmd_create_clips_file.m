function mscmd_create_clips_file(raw_path,clusters_path,output_clips_path,output_index_path,opts)
%MSCMD_CREATE_CLIPS_FILE - Create a couple files that may be used by
%mountainview to view clips
%
%MountainView can be used to view a collection of all clips, organized by
%label (or spike type). However, for memory efficiency, it requires the
%clips file to be prepared with a particular format, along with an index
%file. Compare with ms_create_clips_index.
%
% Syntax:  mscmd_create_clips_file(input_path,clusters_path,output_clips_path,output_index_path,opts)
%
% Inputs:
%    raw_path - path to MxN array of raw or preprocessed data
%    clusters_path - path to file containing times/labels information
%                    according to the particular format outlined in the
%                    docs
%    output_clips_path - path to create the output clips (MxTxNC)
%                        where T=opts.clips_size
%    output_index_path - path to create the index file used by
%                        mountainview. See ms_create_clips_index for more
%                        info.
%
% Other m-files required: mscmd_exe
%
% See also: ms_create_clips_index

% Author: Jeremy Magland
% Jan 2016; Last revision: 13-Feb-2016

if (nargin<5) opts=struct; end;

cmd=sprintf('%s create_clips_file --input=%s --clusters=%s --output=%s --index_out=%s --clip_size=%d ',mscmd_exe,raw_path,clusters_path,output_clips_path,output_index_path,opts.clip_size);

fprintf('\n*** CREATE CLIPS FILE ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end
