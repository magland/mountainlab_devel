function [firingsfile,info] = validspike_wrapper(tsfile,outdir,o)
% VALIDSPIKE_WRAPPER  calls validspike if installed
%
% [firingsfile,info] = validspike_wrapper(tsfile,outdir,o)
%  runs validspike package sorting (greedy glutton fitting, treating all
%  channels together.)
%
% Inputs:
%    tsfile - path to .mda of MxN raw signal timeseries
%    outdir - path to existing DIRECTORY where all output will be written
%    o - (optional) sorting options, see def_sort_opts in this
%                   script
%
% Outputs:
%    firingsfile - path to the firings.mda output file
%    info - struct with fields:
%           filtfile - filtered timeseries
%           prefile - path to the preprocessed timeseries (filt and whitened)
%
% Note: if validspike not found, returns raw timeseries even as filtered/pre
%  and no firings.
%
% See: validspike package.
%
% User must adjust code to point to their compiled validspike installation!

% Barnett 4/11/16

% stuff that's always done to not break things even if validspike absent...
firingsfile = [outdir,'/firings.mda'];

validspikepath = '~/validspike';     % *** USER SET UP THIS ***

if ~exist(validspikepath,'dir')   % exit gracefully with valid but dummy output
  warning('validspike package not installed at known location; check the path in validspike_wrapper!');
  info.prefile = tsfile;         % this is fake
  info.filtfile = tsfile;        % "
  writemda([],firingsfile);
return; end                       % done

savedpath = path;  % set up path
addpath(validspikepath); run('vs_startup');  % also points to spikespy
rmpath([validspikepath,'/../spikespy/matlab'])   % since writemda clash

defo.samplerate = 30000; % default opts
defo.freq_min = 300;   % filter
defo.beta = 3;      % upsampling for det. >3 seems worse for Harris
defo.num_fea = 15;   % PCA
defo.nlps = 10;       % neg log prior for fit
if nargin<3, o=struct; end; o = ms_set_default_opts(o,defo); % setup opts

disp('reading...'); Y = readmda(tsfile);
fprintf('\tM=%d, N=%d (%.3g seconds)\n',size(Y,1), size(Y,2), size(Y,2)/o.samplerate)
disp('filtering...');
d.A = freqfilter(Y,o.samplerate,o.freq_min,[]); % filter and noise-unmix channels
d.samplefreq = o.samplerate; d.dt=1/d.samplefreq;
info.filtfile = [outdir,'/pre0.mda'];
disp('write out filt...'); writemda(d.A,info.filtfile);
disp('whitening...'); d = channelprewhiten(d,[]);
info.prefile = [outdir,'/pre.mda'];
disp('write out filt+whitened...');
writemda(d.A,info.prefile);          % tell output where flit+whitened data is

copts.upsampfac = o.beta;
copts.num_fea=o.num_fea;
sopts.nlps = o.nlps;
[times labels] = spikesort_timeseries(d.A,d.samplefreq,copts,[],sopts); % do it!

disp('write firings...');
writemda([0*times(:)'; times(:)'; labels(:)'],firingsfile);  % dummy channel #s

path(savedpath);  % restore path
