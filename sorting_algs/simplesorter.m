function [firingsfile,prefile]=simplesorter(rawfile,output_dir,o)
% SIMPLERSORTER   Minimal spike sorting code in MATLAB w/ std MS MDA interface.
%
% [firingsfile,prefile]=simplesorter(rawfile,output_dir,opts)
%  runs a simple sorter based on clustering, treating all channels together.
%
% Inputs:
%    rawfile - path to .mda of MxN raw signal data
%    output_dir - path to existing DIRECTORY where all output will be written
%    sort_opts - (optional) sorting options, see def_sort_opts in this
%                   script
%
% Outputs:
%    firingsfile - path to the firings.mda output file
%    prefile - path to the preprocessed raw data file
%
% Other m-files required: isosplit2, processing/ms_*
%
% To test: run without arguments (requires spikespy)

% Barnett based on JFM ds001_sort.m interface 3/16/16. Self-test 3/18/16

if nargin==0, test_simplesorter; return; end

def_sort_opts.clip_size=50;          % all opts defaults
def_sort_opts.samplefreq=30000;
def_sort_opts.freq_min=300;
def_sort_opts.freq_max=10000; %inf;
def_sort_opts.detect_threshold=100;   % in absolute units
def_sort_opts.detect_interval=10;
def_sort_opts.detect_beta=10;
def_sort_opts.num_fea = 10;           % # features

if nargin<3 o=struct; end; o=ms_set_default_opts(o,def_sort_opts); % setup opts

disp('reading...'); Y = readmda(rawfile);
o_filter.samplefreq=o.samplefreq;
o_filter.freq_min=o.freq_min;
o_filter.freq_max=o.freq_max;
disp('filtering...'); Yf = ms_filter(Y,o_filter);
disp('whitening...'); Yf = ms_whiten(Yf);
prefile = [output_dir,'/pre.mda'];
disp('write out filt+whitened...');
writemda(Yf,prefile);          % tell output where flit+whitened data is
o_detect.detect_threshold=o.detect_threshold;
o_detect.detect_interval=o.detect_interval;
o_detect.clip_size=o.clip_size;
o_detect.beta=o.detect_beta;
disp('detect3...'); times = ms_detect3(Y,o_detect);   % note subsample acc
disp('extract clips2...'); clips = ms_extract_clips2(Yf,times,o.clip_size);
disp('features...'); FF = ms_event_features(clips,o.num_fea);
disp('cluster...'); [labels, info] = isosplit2(FF);
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
[firingsfile,~] = simplesorter(d.signal,d.outdir,opts);

% load and view the EC input signal with firings output...
Y = readmda(d.signal);
firings = readmda(firingsfile); times=firings(2,:); labels=firings(3,:);
spikespy({Y,times,labels,'simple sorter'});
K = max(labels); pops = histc(labels,1:K); disp('populations n_k vs k:');
fprintf('\t%d',1:K); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');
