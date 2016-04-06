function d = demo_dataset(forceregen)
% DEMO_DATASET  return struct pointing to demo data files, generate if needed
%
% d = demo_dataset(forceregen) returns struct pointing to a demo timeseries data
%  set.
% Inputs:
%    forceregen: if true, force regeneration w/ random seed (default false)
% Outputs:
%    d struct with fields:
%            d.timeseries - MDA filename for raw EC signal timeseries (signal)
%            d.truefirings - MDA filename for true firings
%            d.truewaveforms - MDA filename for true waveforms
%            d.samplerate - time series sample rate in samples/sec
%            d.outdir - path to a working directory to be used for output files
%
% The only self-test is to run without arguments
%
% See also: LOADDEMOWAVEFORMS

% Barnett 3/18/16
mfile_path=fileparts(mfilename('fullpath'));
demodir = [mfile_path,'/demo_data'];
if ~exist(demodir,'dir'), mkdir(demodir); end
d.outdir = [demodir,'/output'];
if ~exist(d.outdir,'dir'), mkdir(d.outdir); end

d.timeseries = [demodir,'/demotimeseries.mda'];
d.truefirings = [demodir,'/demotruefirings.mda'];
d.truewaveforms = [demodir,'/demotruewaveforms.mda'];
d.name = 'demo data from EJ 2005-04-26 elec359 K=7';
d.samplerate = 2e4;           % samples/sec to generate

if nargin<1, forceregen=0; end
if forceregen | ~exist(d.timeseries,'file')  % could check other out files too?
  disp('generating demo timeseries...')
  info = gen_demo_data(d);
end
%%%%%%%%%%%
  
function info = gen_demo_data(d)
% GEN_DEMO_DATA  generate and write out demo data, groundtruthed
%
% Params baked in. Specific to EJ K7 data.  input is a struct

% Barnett 3/18/16, simplified from write_demo_timeseries_func

[W samplerateW] = loaddemowaveforms;         % should switch different W's
writemda32(W,d.truewaveforms);           % in case the validator needs them
o_synth.upsamplefac = round(samplerateW/d.samplerate); % allows upsampled W
K = size(W,3);
T = 120;         % desired length in seconds
N = round(T*d.samplerate);         % number of time points
rates = 0.3*2.^(0:K-1);     % range of firing rates, gives up to >10 Hz
rates = rates/d.samplerate;   % convert to per sample
o_firings.amplsig = 0.0;  % 0.2       % relative amplitude std dev about 1
[times labels ampls] = randomfirings(N,rates,o_firings);
peakchans = 0*times;
% write out ground-truth firings in correct format... since validator needs
writemda64([peakchans;times;labels;ampls],d.truefirings);
tic; Y = synthesize_timeseries(W,N,times,labels,ampls,o_synth); toc
eta = 20;               % noise std deviation per sample per channel
Y = Y + eta * randn(size(Y));
writemda32(Y,d.timeseries);
info = [];    % dummy for now