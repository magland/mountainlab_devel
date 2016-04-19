function d = grab_martinez2009_dataset(n)
% Get Martinez-Quiroga 2009 JNM 1-chan sim datasets w/ known firings.
%
% d = grab_martinez2009_dataset(n) where n = 1...5 returns struct pointing
% to MDA conversion from one of 5 datasets.

% Barnett 4/11/16

if nargin<1, n=1; end   % default

d.outdir = '/tmp/output'; if ~exist(d.outdir,'dir'), mkdir(d.outdir); end  % fix

mfile_path=fileparts(mfilename('fullpath'));
ext = [mfile_path,'/../ext_datasets'];
dir = [ext,'/MartinezQuiroga2009_sims'];
rawfile = sprintf('simulation-%d',n);
fname = strcat(dir,'/',rawfile);
fprintf('reading %s ...\n',fname)
load(fname);
d.samplerate = 1e3/samplingInterval;                  % was in ms
d.timeseries = [d.outdir,'/martinez2009_raw.mda'];
writemda(data, d.timeseries,'float32');
d.name = sprintf('Martinez09-sim%d',n);
truelabels = 1+spike_class{1};        % convert from 0-indexed labels
K = max(truelabels);
fprintf('dataset has K=%d true labels\n',K)
timeoff = 22.0;    % seems like firing times were # samples behind actual peaks
truetimes = spike_times{1} + timeoff;
truefirings = [0*truetimes; truetimes; truelabels];  % chan; time; label
d.truefirings = [d.outdir,'/martinez2009_truefirings.mda'];
writemda(truefirings, d.truefirings,'float64');
