function [clips,templates] = reorderbytemplatenorm(firingsfile,signalfile,o)
% helper for postprocessing. todo: docs
%
% Note: This file rewrites the firingsfile!
%
% Outputs:
%   clips - in case you want them
%   templates - in the new sorted order

% Barnett 3/21/16 inspired by validspike/stageB/spikesort_clips.m

if nargin<3, o=[]; end
if ~isfield(o,'clip_size'), o.clip_size = 60; end

firings = readmda(firingsfile); times=firings(2,:); labels=firings(3,:);
Y = readmda(signalfile);
clips = ms_extract_clips2(Y,times,o.clip_size);  % get raw clips & wf's
templates = ms_templates(clips,labels);
[M T K] = size(templates);
Wnrms = sqrt(sum(reshape(templates.^2,[M*T K]),1));
[~,perm] = sort(Wnrms,'descend');
[~,iperm] = sort(perm);                 % inverse of permutation
templates = templates(:,:,perm);        % apply the perm to meaningful labels...
ii=(labels>0 & labels<=K); labels(ii) = iperm(labels(ii));
firings(3,:) = labels;
writemda(firings,firingsfile,'float64');
