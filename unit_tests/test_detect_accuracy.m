function test_detect_accuracy(detfunc,o_det,forceregen)
% TEST_DETECT_ACCURACY  measure accuracy of any ms_detect function
%
% test_detect_accuracy(detfunc,o_det,forceregen) tests any detection function
%  on synthetic data that has known firing times. This is useful for assessing
%  ability to measure peak time to subsample accuracy.
%
% Inputs:
%  detfunc - handle to detection function
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

function demo_test_detect_accuracy
o.detect_threshold = 90;   % absolute (uV) units, for EJ data
o.detect_interval = 5;
o.clip_size = 30;          % only affects ends of timeseries
regendata = 1;             % toggle this as you please
% ahb trying various detection algs...
fprintf('\nold detect:\n')             % kept in scratch_ahb
v = path; addpath scratch_ahb
test_detect_accuracy(@ms_detect,o,regendata); title('old detect');
fprintf('\ndetect3:\n')
test_detect_accuracy(@ms_detect3,o); title('detect3');
o.jiggle = 1; fprintf('\ndetect4 jiggle=1:\n')
test_detect_accuracy(@ms_detect4,o); title('detect4 jiggle=1');
o.jiggle = 2; fprintf('\ndetect4 jiggle=2:\n')
test_detect_accuracy(@ms_detect4,o); title('detect4 jiggle=2');
o.num_features=10; fprintf('\ndetect4 jiggle=2 numfea=10:\n')
test_detect_accuracy(@ms_detect4,o); title('detect4 jiggle=2 numfea=30');
path(v);
% seems like jiggle=1 helps, but not any higher, and numfea around 15 best
% Also, adding jiggle to single times set is slightly worse than appending.
