function d = grab_neurocube_dataset(n)
% gets Camunas-Mesa-Quiroga '14 NeuroCube synthetic tetrode w/ known firings
%
% d = grab_neurocube_dataset(n)
%  n=1:   30 sec
%  n=2:   300 sec = 5 min
% Barnett 4/26/16

if nargin==0, n=1; end

% keep output dir in a fast-IO drive:
d.outdir = [tempdir,'output']; if ~exist(d.outdir,'dir'), mkdir(d.outdir); end

mfile_path=fileparts(mfilename('fullpath'));
ext = [mfile_path,'/../ext_datasets'];
dir = [ext,'/NeuroCube/mat'];

d.samplerate = 24000;        % see Par_sim.sr in neurocube.m
d.timeseries = [d.outdir,'/neurocube_default_tet.mda'];
d.truefirings = [d.outdir,'/neurocube_default_tet_truefirings.mda'];

if n==1
  load([dir,'/Spikesim_default_tet_4-25-16.mat']);
  d.name = 'NeuroCube-default-tet';
elseif n==2
  load([dir,'/Spikesim_tet_5min.mat']);
  d.name = 'NeuroCube-tet-5min';
end
writemda32(data,d.timeseries);
Wthresh = 20.0;                   % template abs threshold
gnd = find(max(abs(Close_Spikeshapes),[],2)>Wthresh); % labels to call "ground truth"
% a few neurons that have large peaks
times = []; labels = [];
% relabel them 1,..,K and output as ground truth...
for k=1:numel(gnd)
  ev = find(Close_neurons(:,1)==gnd(k));  % events of this type
  times = [times Close_neurons(ev,2)'];   % row append
  labels = [labels k*ones(1,numel(ev))];  % "
end
times = d.samplerate*times/1e3;    % convert ms -> samples
times = times + 21.0;              % apparent offset of 21 samples (as in 2009!)
writemda64([0*times; times; labels],d.truefirings);  % peakchans=0


if 0 % plotting to see if any correlation of Y clips with firings:
  spikespy({data,times,labels,'NeuroCube'});
  K = numel(gnd);
  pops = histc(labels,1:K); [1:K; pops]
  figure; ms_view_templates(ms_templates(ms_extract_clips2(data,times,60),labels),struct('pops',pops));
end