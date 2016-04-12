function [firingsfile,info]=simplesorter(rawfile,output_dir,o)
% SIMPLERSORTER   Minimal spike sorting code in MATLAB w/ std MS MDA interface.
%
% [firingsfile,info] = simplesorter(rawfile,output_dir,o)
%  runs a simple sorter based on clustering, treating all channels together.
%
% Inputs:
%    rawfile - path to .mda of MxN raw signal data
%    output_dir - path to existing DIRECTORY where all output will be written
%    o - (optional) sorting options, see def_sort_opts in this
%                   script
%
% Outputs:
%    firingsfile - path to the firings.mda output file
%    info - struct with fields:
%           filtfile - filtered timeseries
%           prefile - path to the preprocessed timeseries (filt and whitened)
%
% Other m-files required: isosplit2, processing/ms_*
%
% To test: run without arguments (requires spikespy)

% Barnett based on JFM ds001_sort.m interface 3/16/16. Self-test 3/18/16
% added polarity 3/21/16, added detect4 4/7/16

if nargin==0, test_simplesorter; return; end

def_sort_opts.clip_size=50;          % all opts defaults
def_sort_opts.samplerate=30000;
def_sort_opts.freq_min=300;
def_sort_opts.freq_max=10000; %inf;

def_sort_opts.detect_threshold=3.5;   % in std dev units
def_sort_opts.detect_interval=10;
def_sort_opts.detect_beta=10;
def_sort_opts.detect_meth='p';
def_sort_opts.detect_polarity='b';

def_sort_opts.num_fea = 10;           % # features
def_sort_opts.verb = 0;                % verbosity

if nargin<3 o=struct; end; o=ms_set_default_opts(o,def_sort_opts); % setup opts

disp('reading...'); Y = readmda(rawfile);
fprintf('\tM=%d, N=%d (%.3g seconds)\n',size(Y,1), size(Y,2), size(Y,2)/o.samplerate)
o_filter.samplerate=o.samplerate;
o_filter.freq_min=o.freq_min;
o_filter.freq_max=o.freq_max;
disp('filtering...'); Yf = ms_bandpass_filter(Y,o_filter);
info.filtfile = [output_dir,'/pre0.mda'];
disp('write out filt...'); writemda(Yf,info.filtfile);
disp('whitening...'); Yf = ms_whiten(Yf);   % normalizes
info.prefile = [output_dir,'/pre.mda'];
disp('write out filt+whitened...');
writemda(Yf,info.prefile);          % tell output where flit+whitened data is
o_detect.detect_threshold=o.detect_threshold;
o_detect.detect_interval=o.detect_interval;
o_detect.clip_size=o.clip_size;
o_detect.beta=o.detect_beta;
o_detect.meth = o.detect_meth;
o_detect.polarity = o.detect_polarity;
disp('detect4...'); times = ms_detect4(Yf,o_detect);   % note subsample acc
fprintf('\tfound %d events\n',numel(times))
disp('extract clips2...'); clips = ms_extract_clips2(Yf,times,o.clip_size);
disp('features...'); FF = ms_event_features(clips,o.num_fea);
disp('cluster...'); [labels, info.isosplit2info] = isosplit2(FF);
if o.verb, figure; ms_view_clusters(FF,labels); end
firingsfile = [output_dir,'/firings.mda'];
disp('write firings...');
writemda([0*times; times; labels],firingsfile);   % use dummy channel #s
K = max(labels);

%pops = histc(labels,1:K);
%fprintf('populations (sorter ordering):\n'); fprintf('\t%d',pops); fprintf('\n');

%%%%%%%%
function test_simplesorter  % taken from examples/driver_simplesorter

d = demo_dataset;   % get struct pointing to demo data files
opts.samplerate = d.samplerate;
[firingsfile,~] = simplesorter(d.timeseries,d.outdir,opts);

% load and view the EC input signal with firings output...
Y = readmda(d.timeseries);
firings = readmda(firingsfile); times=firings(2,:); labels=firings(3,:);
if exist('spikespy')
  spikespy({Y,times,labels,'simple sorter'});
end
K = max(labels); pops = histc(labels,1:K); disp('populations n_k vs k:');
fprintf('\t%d',1:K); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');
