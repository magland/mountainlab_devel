function test_detect_accuracy(detfunc,o_det,forceregen)
% TEST_DETECT_ACCURACY  accuracy of any ms_detect function
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

% Barnett 3/3/16, 3/11/16

if nargin==0, demo_test_detect_accuracy; return; end
if nargin<2, o_det = []; end
if nargin<3, forceregen = 0; end

[Yfile truefiringsfile trueWfile o.samplefreq] = get_default_dataset([],forceregen); % make demo data and return names of its files

mfile_path=fileparts(mfilename('fullpath'));
outdir=[mfile_path,'/demo_data/output'];
if ~exist(outdir,'dir') mkdir(outdir); end;

Y=readmda(Yfile);                 % use raw

%o.freq_min = 100; o.freq_max = 5000; Y = ms_filter(Y,o); % doesn't help

times = detfunc(Y,o_det);   % do the thing

% compare to true
truefirings=readmda(truefiringsfile);       % read sort output files
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
o.detect_threshold = 90;
o.detect_interval = 5;
o.clip_size = 30;          % only affects ends of timeseries
test_detect_accuracy(@ms_detect,o); title('old detect'); % ,1) if want regen data
test_detect_accuracy(@ms_detect3,o); title('new detect3');
