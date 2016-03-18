% Example script showing how to synthesize random time series and write to MDA.

% Barnett 2/19/16; allows upsampled waveforms 2/25/16. duplicates demo_dataset

[W samplerateW] = loaddemowaveforms;
samplerate = 2e4;           % in samples/sec
o_synth.upsamplefac = round(samplerateW/samplerate); % allows upsampled W
K = size(W,3);
T = 120;                    % desired length in seconds
N = round(T*samplerate);    % number of time points
rates = 0.3*2.^(0:K-1);     % range of firing rates, gives up to >10 Hz for K=7
rates = rates/samplerate;   % convert to per sample
o_firings.amplsig = 0.2;    % relative amplitude std dev about 1
[times labels ampls] = randomfirings(N,rates,o_firings);
peakchans = 0*times;
% write out ground-truth firings in correct format...
writemda([peakchans;times;labels;ampls],'unit_tests/demo_data/truefirings.mda');
Y = synthesize_timeseries(W,N,times,labels,ampls,o_synth);
eta = 20;               % noise std deviation per sample per channel
Y = Y + eta * randn(size(Y));
writemda(Y,'unit_tests/demo_data/exampletimeseries.mda');

% take a look...
%spikespy({Y,round(times),labels,'demo Y'});
%figure; tmax=1e4; plot(Y(:,1:tmax)','.-'); hold on; vline(times(times<tmax));
