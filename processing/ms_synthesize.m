function Y = ms_synthesize(W,N,times,labels,ampls,opts)
% MS_SYNTHESIZE.  Make timeseries via linear addition of spikes forward model
%
% Y = ms_synthesize(W,T,times,labels,ampls,opts) outputs a synthetic time series
% given waveforms and firing information; ie, it applies the forward model.
% No noise is added.
%
% Inputs:
%  W - waveform (templates) array, M by T by K. (note: assumed not upsampled
%      unless opts.upsamplefac > 1).
%  N - integer number of time points to generate.
%  times - (1xL) list of firing times, integers or real,
%          where 1st output time is 1 and last is N (ie, sample index units).
%  labels - (1xL) list of integer labels in range 1 to K.
%  ampls - (optional, 1xL) list of firing amplitudes (if absent or empty,
%          taken as 1).
%  opts - (optional struct) controls various options such as
%         opts.upsamplefac : integer ratio between sampling rate of W and for
%                            the output timeseries.
% Outputs:
%  Y - (MxN, real-valued) timeseries

% Barnett 2/19/16 based on validspike/synthesis/spikemodel.m
% todo: faster C executable acting on MDAs I/O.

if nargin==0, test_ms_synthesize; return; end
if nargin<5 || isempty(ampls), ampls = 1.0+0*times; end      % default ampl
if nargin<6, opts = []; end
if ~isfield(opts,'upsamplefac'), opts.upsamplefac = 1; end

[M T K] = size(W);
tcen = floor((T+1)/2);    % center firing time in waveform
L = numel(times);
if numel(labels)~=L, error('times and labels must have same # elements!'); end
Y = zeros(M,N);
for j=1:L            % loop over spikes adding in each one
  toff = round(times(j)) - tcen;   % integer offset of start of W from output Y
  i = max(1,1-toff):min(T,N-toff);      % valid indices in waveform's 1:T
  Y(:,i+toff) = Y(:,i+toff) + ampls(j)*W(:,i,labels(j));
end

function test_ms_synthesize
% make simple variable-width Gaussian waveforms...
M = 4;   % # channels
T = 30;  % # timepoints for waveform
K = 5;   % # neuron types
W = zeros(M,T,K);
tcen = floor((T+1)/2); % center index
t = (1:T) - tcen;          % offset time grid
for k=1:K
  wid = 3*exp(0.5*randn);         % Gaussian width in time
  pulse = exp(-0.5*t.^2/wid^2);
  W(:,:,k) = randn(M,1) * pulse;            % outer prod
end
% make firing info...
N = 1e6;   % total time points
L = 1e4;   % number of firings
times = round(rand(1,L)*(N-1)+1);  % integer for now
labels = randi(K,1,L);
tic
Y = ms_synthesize(W,N,times,labels);
toc
spikespy({Y,times,labels,'test ms_synthesize'});
