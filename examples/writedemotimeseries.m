% example script showing how to synthesize random time series and write to MDA.
% Note: this script is used by examples/sort_demo.m
% Barnett 2/19/16

clear
[W samplerate] = loaddemowaveforms;
K = size(W,3);
T = 120;         % desired length in seconds
N = round(T*samplerate);         % number of time points
% mean firing rates in Hz, K=7 must match demo waveforms:
rates = 0.1*2.^(0:K-1);     % range of firing rates
rates = rates/samplerate;   % convert to per sample
o_firings.amplsig = 0.2;         % amplitude variation
[times labels ampls] = randomfirings(N,rates,o_firings);
peakchans = 0*times;
% write out ground-truth firings in correct format...
writemda([peakchans;times;labels;ampls],'unit_tests/demo_data/truefirings.mda');
Y = ms_synthesize(W,N,times,labels,ampls);
eta = 20;               % noise std deviation per sample per channel
Y = Y + eta * randn(size(Y));
writemda(Y,'unit_tests/demo_data/demotimeseries.mda');
%spikespy({Y,times,labels,'demo Y'});  % if you want a look
