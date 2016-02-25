function raw_path=default_test_dataset(name)

if nargin<1 name='synthetic1'; end;

mfile_path=fileparts(mfilename('fullpath'));

if (strcmp(name,'synthetic1'))
    raw_path=[mfile_path,'/demo_data/synthetic1.mda'];
    write_demo_timeseries(raw_path);
elseif (strcmp(name,'synthetic2'))
    raw_path=[mfile_path,'/demo_data/synthetic1.mda'];
    write_demo_timeseries(raw_path);
else
    error('Unknown name: %s',name);
end;

end



function write_demo_timeseries(path)
% Example script showing how to synthesize random time series and write to MDA.
% Note: this script is used by examples/sort_demo.m

% Barnett 2/19/16; allows upsampled waveforms 2/25/16

[W samplerateW] = loaddemowaveforms;
samplerate = 2e4;           % in samples/sec
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
% write out ground-truth firings in correct format...
%writemda([peakchans;times;labels;ampls],'unit_tests/demo_data/truefirings.mda');
tic; Y = ms_synthesize(W,N,times,labels,ampls,o_synth); toc
eta = 20;               % noise std deviation per sample per channel
Y = Y + eta * randn(size(Y));
writemda_if_necessary(Y,path);

end
