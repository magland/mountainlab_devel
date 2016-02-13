function clips=ms_extract_clips(Y,times,clip_size)
%MS_EXTRACT_CLIPS - Extract clips from raw/preprocessed data array centered
%                   at a collection of event times.
%
%Consider using mscmd_extract_clips
%
% Syntax:  [clips] = ms_extract_clips(Y,times,clip_size)
%
% Inputs:
%    Y - MxN array of raw or pre-processed data
%    times - 1xNC array of integer timepoints (event times)
%    clip_size - integer, e.g. 100 (timepoints)
%
% Outputs:
%    clips - MxTxNC array of extracted clips
%
% Example:
%    clips=ms_extract_clips(X,times0,100);
%    FF=ms_event_features(clips,3);
%    labels=isosplit2(FF)
%    ms_view_clusters(FF,labels);
%
% Other m-files required: none
%
% See also: mscmd_extract_clips, ms_event_features, spikespy

% Author: Jeremy Magland
% Jan 2016; Last revision: 13-Feb-2016

[M,N]=size(Y);
T=clip_size;
C=length(times);

clips=zeros(M,T,C);
tt1=-ceil((clip_size)/2);
tt2=tt1+clip_size-1;
inds=find((times+tt1>=1)&(times+tt2<=N));
%if (min(times+tt1)<1) error('Invalid time in extract_clips'); end;
%if (max(times+tt2)>N) error('Invalid time in extract_clips'); end;
%for j=1:C
for j=inds
	clips(:,:,j)=Y(:,times(j)+tt1:times(j)+tt2);
end;

end
