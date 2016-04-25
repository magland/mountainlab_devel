function d = grab_testmerge_dataset()
% struct pointing to simple equi-peak synth dataset useful for testing merge stage.
%
% This writes a synth dataset with 2 channels on which a single neuron type
%  fires with equal peaks, but with a time shift. Useful to test merging across
%  channels. No upsampling. Always regenerated, random Poisson firings.
%
% See also: ms_merge_across_channels.    Barnett 4/25/16

d.outdir = [tempdir,'/testmerge'];
if ~exist(d.outdir,'dir'), mkdir(d.outdir); end
d.name = 'testmerge-equalpeaks';

dt = 5;   % now draw from ms_merge_across_channels self-test to make waveforms:
T = 30;
M=2;
K=1;           % one type with peak similar onb
templates = nan(M,T,K);
Tcen = floor((T+1)/2);  % defn of center time of clip
w0=3;   % Gaussian width in samples (w0=3 makes r<0.5 for dt>=5)
t = (1:T)-Tcen;
templates(:,:,1) = [1.0*exp(-.5*t.^2/w0^2); 1.0*exp(-.5*(t-dt).^2/w0^2)];
N = 1e7;   % total duration
L = 1e4;   % roughly how many true firings
truetimes = randomfirings(N,L/N,struct('refractory',30));  % have refractory period
% (it's interesting how missing this period messes stuff up).

d.truefirings = [d.outdir,'/truefirings.mda'];
truelabels = 1+0*truetimes;   % all labels 1
writemda64([0*truetimes;truetimes;truelabels],d.truefirings);
Y = synthesize_timeseries(templates,N,truetimes,truelabels);
eta = 0.25;               % noise std deviation per sample per channel
Y = Y + eta * randn(size(Y));
d.timeseries = [d.outdir, '/raw.mda'];
writemda32(Y,d.timeseries);
d.samplerate = 2e4;          % this is a dummy, since testing stuff indep of this
