function [times labels ampls] = randomfirings(N,rates,opts)
% RANDOMFIRINGS.  Make Poisson random firing times and labels.
%
% [times labels ampls] = randomfirings(N,rates,opts) makes independent random
%   Poissonian firing trains at the requested mean rates.
%
% Inputs:
%   N - maximum time (integer sample units)
%   rates (1xK) - mean firing rates for each neuron type, in per sample units
%   opts - (optional) struct controlling options:
%          opts.amplsig = std deviation of amplitudes about 1.0 (default 0)
% Outputs:
%   times - (1xL) real-valued firing times in range 1 to N
%   labels - (1xL) labels in 1...K
%   ampls - (1xL) firing amplitudes (reals).
%
% simplified from validspike synth_Poissonspiketrain.   Barnett 2/19/16

if nargin==0, test_randomfirings; return; end
if nargin<3, opts=[]; end
if ~isfield(opts,'amplsig'), opts.amplsig = 0; end

K = numel(rates);
times = []; labels = []; ampls = [];
tmin = 1; tmax = N;                      % allowed time range
for k=1:K
  tk = []; t = tmin - log(rand(1))/rates(k);  % generate first firing
  while (t<tmax)
    tk = [tk t]; t = t - log(rand(1))/rates(k); % subsequent firings
  end
  Ne(k) = numel(tk);
  times = [times tk];
  labels = [labels k*ones(1,Ne(k))];
  ampls = [ampls, 1.0 + opts.amplsig*randn(1,Ne(k))];
end

function test_randomfirings
opts.amplsig = 0.2;
[times labels ampls] = randomfirings(1e4,[1e-1 1e-2 1e-3],opts);
fprintf('populations: '); histc(labels,1:3)
hist(ampls,50); title('firing ampls');
figure; showfirings(times,labels); title('firing times');
spikespy({times,labels});
