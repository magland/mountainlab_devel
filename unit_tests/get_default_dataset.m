function [Yfile truefiringsfile trueWfile samplerate] = get_default_dataset(name,forceregen)
% GET_DEFAULT_DATASET.  Generate synthetic dataset and return path
%
% [Yfile truefiringsfile trueWfile samplerate] = default_test_dataset(name,forceregen)
%
% This is the universal function to call to get a synthetic datasets.
%
% Input:
%  name - a string choosing the test dataset type (if absent, default used)
%  forceregen (optional) - flag: if 0 (default) only regenerate data files if
%                                     absent. if 1, force regeneration.
%
% Outputs: first 3 are filenames generated of Y timeseries, and truth firings and
%  waveforms. 4th output is sample rate in samples/sec.

% Magland+Barnett 2/25/16. forceregen 3/11/16

if nargin<1 || isempty(name) || strcmp(name,''), name='EJ K7'; end
if nargin<2, forceregen=0; end

mfile_path=fileparts(mfilename('fullpath'));

if (strcmp(name,'EJ K7'))
  head=[mfile_path,'/demo_data/EJK7'];  % relative to this function!
  [Yfile truefiringsfile trueWfile samplerate] = write_demo_timeseries_func(head,forceregen);
elseif (strcmp(name,'synthetic2'))
  % etc ...
else
  error('Unknown name: %s',name);
end

%%%%%%%%%%%

function [Yfile truefiringsfile trueWfile samplerate] = write_demo_timeseries_func(head,forceregen)
% Function packaging the synthesis script for the "demo" waveforms.
%
% Example script showing how to synthesize random time series and write to MDA.

% Barnett 2/19/16; allows upsampled waveforms 2/25/16. forceregen flag 3/11/16

Yfile = [head '.mda'];
truefiringsfile = [head '_truefirings.mda'];
trueWfile = [head '_trueW.mda'];
samplerate = 2e4;           % in samples/sec

if forceregen || ~exist(Yfile)    % lack of Yfile is enough to trigger regen...
  [W samplerateW] = loaddemowaveforms;         % should switch different W's
  writemda32(W,trueWfile);           % in case the validator needs them
  o_synth.upsamplefac = round(samplerateW/samplerate); % allows upsampled W
  K = size(W,3);
  T = 120;         % desired length in seconds
  N = round(T*samplerate);         % number of time points
  % mean firing rates in Hz, K=7 must match demo waveforms:
  rates = 0.3*2.^(0:K-1);     % range of firing rates, gives up to >10 Hz
  rates = rates/samplerate;   % convert to per sample
  o_firings.amplsig = 0.0;  % 0.2       % amplitude variation
  [times labels ampls] = randomfirings(N,rates,o_firings);
  peakchans = 0*times;
  % write out ground-truth firings in correct format... since validator needs
  writemda64([peakchans;times;labels;ampls],truefiringsfile);
  tic; Y = ms_synthesize(W,N,times,labels,ampls,o_synth); toc
  eta = 20;               % noise std deviation per sample per channel
  Y = Y + eta * randn(size(Y));
  writemda32(Y,Yfile);
end
