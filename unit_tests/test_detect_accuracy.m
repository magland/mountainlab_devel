function test_detect_accuracy(detfunc,o_det,forceregen)
% TEST_DETECT_ACCURACY  measure accuracy of any ms_detect function
%
% test_detect_accuracy(detfunc,o_det,forceregen) tests any detection function
%  on synthetic data that has known firing times. This is useful for assessing
%  ability to measure peak time to subsample accuracy.
%
% Inputs:
%  detfunc - handle to detection function w/ interface: times = detfunc(Y,o_det)
%  o_det - options struct for detection function
%  forceregen (optional) - flag stating whether to force regeneration of demo
%        data (default 0).
%
% Example usage: see enclosed code demo_test_detect_accuracy
%
% Also see: MS_DETECT

% todo: make standard output values (MAE, frac of correct found, ...)

% Barnett 3/3/16, 3/11/16, used demo_dataset & jiggle 3/25/16

if nargin==0, demo_test_detect_accuracy; return; end
if nargin<2, o_det = []; end
if nargin<3, forceregen = 0; end

d = demo_dataset(forceregen);
Y=readmda(d.timeseries);                 % use raw

%o.freq_min = 100; o.freq_max = 5000; Y = ms_filter(Y,o); % doesn't help

times = detfunc(Y,o_det);   % do the thing

% compare to true
truefirings=readmda(d.truefirings);       % read sort output files
truetimes=truefirings(2,:); truelabels=truefirings(3,:);
maxdt = 5;
CC = ms_cross_correlograms([truetimes, times], [1+0*truetimes, 2+0*times],maxdt);
tdiffs = CC{1,2};
figure; hist(tdiffs,100); axis([-maxdt,maxdt,ylim]);
fprintf('mean abs detection time error in range [-maxdt,maxdt]: %.3g samples\n', mean(abs(tdiffs)))
fprintf('rms detection time error in range [-maxdt,maxdt]:      %.3g samples\n\n', sqrt(mean(tdiffs.^2)))

o.max_matching_offset = 2;    % in samples
o.verb = 0;
[~,Q,acc,t] = times_labels_accuracy(truetimes,1+0*truelabels,times,1+0*times,o);
%disp('extended confusion matrix (1st col = spike vs 2nd col = no-spike):'); Q

%%%%%%%%%%%%%%%%%%% 

%%% helpers for the demo...

function fname = fnameify32(X,outdir)
% FNAMEIFY  if array, writes to file and returns filename, otherwise keeps name

% v crude for now.

if nargin<2, outdir=tempdir; end
if ischar(X) || isstring(X)
  fname = X;
else
  fname = [outdir,'/',num2str(randi(1e10)),'.mda'];  % random filename
  writemda32(X,fname);
end

function times = run_mscmd_detect(Y,o);
outdir = tempdir; detectfile = [outdir,'/detect.mda']; 
mscmd_detect(pathify32(Y),detectfile,o);
firings = readmda(detectfile);
times = firings(2,:);

function times = run_mscmd_detect3(Y,o);
outdir = tempdir; detectfile = [outdir,'/detect.mda']; 
mscmd_detect3(pathify32(Y),detectfile,o);               % *** broken now
firings = readmda(detectfile);
times = firings(2,:);

%%%%%%%

function demo_test_detect_accuracy
o.detect_threshold = 90;   % absolute (uV) units, for EJ data
o.detect_interval = 5;
o.clip_size = 30;          % only affects ends of timeseries
o.polarity='m';
regendata = 1;             % toggle this as you please
% ahb trying various detection algs...
fprintf('\nold ms detect:\n')             % kept in scratch_ahb
v = path; addpath scratch_ahb
test_detect_accuracy(@ms_detect,o,regendata); title('old ms detect'); drawnow
fprintf('\nms detect3:\n')
test_detect_accuracy(@ms_detect3,o); title('ms detect3'); drawnow
o.jiggle = 1; fprintf('\nms detect4 jiggle=1:\n')
test_detect_accuracy(@ms_detect4,o); title('ms detect4 jiggle=1'); drawnow
o.jiggle = -1; fprintf('\nms detect4 jiggle=-1 (increase jiggle):\n')
test_detect_accuracy(@ms_detect4,o); title('ms detect4 jiggle=-1'); drawnow
o.jiggle = 2; fprintf('\nms detect4 jiggle=2:\n')
test_detect_accuracy(@ms_detect4,o); title('ms detect4 jiggle=2'); drawnow
o.num_features=10; fprintf('\nms detect4 jiggle=2 numfea=10:\n')
test_detect_accuracy(@ms_detect4,o); title('ms detect4 jiggle=2 numfea=10'); drawnow
o.sign=-1; o.individual_channels = 0;
test_detect_accuracy(@run_mscmd_detect,o); title('mscmd detect'); drawnow
o.sign=-1; o.upsampling_factor = 10; o.num_pca_denoise_components = 5;
o.pca_denoise_jiggle = 2;
o.individual_channels = 0;   % jfm to implement (1 find duplicates of spikes!)
test_detect_accuracy(@run_mscmd_detect3,o); title('mscmd detect3'); drawnow
path(v);
% seems like jiggle=1 helps, but not any higher, and numfea around 15 best
% Also, adding jiggle to single times set is slightly worse than appending.
