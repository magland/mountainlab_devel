function test_detect_accuracy(detfunc,o_det)
% TEST_DETECT_ACCURACY  accuracy of any ms_detect function
%
% test_detect_accuracy(detfunc,o_det)
%
% Example usage: see enclosed code demo_test_detect_accuracy
%
% Also see: MS_DETECT

% Barnett 3/3/16

if nargin==0, demo_test_detect_accuracy; return; end
if nargin<2, o_det = []; end
[Yfile truefiringsfile trueWfile o.samplefreq] = get_default_dataset; % demo data

mfile_path=fileparts(mfilename('fullpath'));
outdir=[mfile_path,'/demo_data/output'];
if ~exist(outdir,'dir') mkdir(outdir); end;

Y=readmda(Yfile);                 % use raw

%o.freq_min = 100; o.freq_max = 5000; Y = ms_filter(Y,o); % doesn't help

times = ms_detect(Y,o_det);

% compare to true
truefirings=readmda(truefiringsfile);       % read sort output files
truetimes=truefirings(2,:); truelabels=truefirings(3,:);
maxdt = 5;
CC = ms_cross_correlograms([truetimes, times], [1+0*truetimes, 2+0*times],maxdt);
tdiffs = CC{1,2};
figure; hist(tdiffs,100);
fprintf('mean abs detection time error in range [-maxdt,maxdt]: %.3g samples\n', mean(abs(tdiffs)))
fprintf('rms detection time error in range [-maxdt,maxdt]:      %.3g samples\n\n', sqrt(mean(tdiffs.^2)))

o.max_matching_offset = 2;    % in samples
[~,Q,acc,t] = times_labels_accuracy(truetimes,1+0*truelabels,times,1+0*times,o);

%%%%%%%%%%%%%%%%%%% 


function demo_test_detect_accuracy
o.detect_threshold = 80;
o.detect_interval = 20;
o.clip_size = 50;
test_detect_accuracy(@ms_detect,o)
test_detect_accuracy(@ms_detect3,o)
