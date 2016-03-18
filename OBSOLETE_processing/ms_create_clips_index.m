function [clips_out,clips_index]=ms_create_clips_index(clips,labels)
%CREATE_CLIPS_INDEX - Create a couple arrays that may be used by mountainview
%
%MountainView can be used to view a collection of all clips, organized by
%label (or spike type). However, for memory efficiency, it requires the
%clips file to be prepared with a particular format, along with an index
%file.
%
% Syntax:  [clips_out,clips_index] = ms_create_clips_index(clips,labels)
%
% Inputs:
%    clips - MxTxNC array of clips (time windows from raw/preprocessed data)
%    labels - 1xNC vector of integer labels
%
% Outputs:
%    clips_out - MxTxNC array of clips organized by label
%    clips_index - 1x(K+1) index vector, used by MountainView
%                  where K=max(labels)
%
% Example: 
%    [clips0,index0]=create_clips_index(clips,labels);
%    writemda('clips0.mda',clips0);
%    writemda('index0.mda',index0);
%    % Now run MountainView and point to these files
%
% Other m-files required: none
%
% See also: ms_extract_clips, ms_mountainview, ms_view_templates_from_clips

% Author: Jeremy Magland
% Jan 2016; Last revision: 13-Feb-2016

[M,T,NC]=size(clips);
clips_out=zeros(M,T,NC);
K=max(labels);
clips_index=zeros(1,K+1);
ind0=0;
for k=1:K
    inds_k=find(labels==k);
    clips_out(:,:,(ind0+1):(ind0+length(inds_k)))=clips(:,:,inds_k);
    clips_index(k)=ind0;
    ind0=ind0+length(inds_k);
end;
clips_index(end)=NC;

end
