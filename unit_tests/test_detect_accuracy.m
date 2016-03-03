function test_detect_accuracy(detfunc,o_det)
% TEST_DETECT_ACCURACY  accuracy of any ms_detect function
%
% Example usage:
%  o.detect_threshold = 100;
%  o.detect_interval=20;
%  o.clip_size=50;
%  test_detect_accuracy(@ms_detect,o)
%  where o_detect is (optional) struct passed to the detector function.
%
%  For now, just reports # misses and false pos, using fixed matching offset.
%
% To do: make report sub-sample t_j error histogram
%
% Also see: MS_DETECT

% Barnett 3/3/16

if nargin<2, o_det = []; end
[Yfile truefiringsfile trueWfile samplerate] = get_default_dataset; % demo data

mfile_path=fileparts(mfilename('fullpath'));
outdir=[mfile_path,'/demo_data/output'];
if ~exist(outdir,'dir') mkdir(outdir); end;

Y=readmda(Yfile);                 % use raw
times = ms_detect(Y,o_det);

% compare to true
truefirings=readmda(truefiringsfile);       % read sort output files
truetimes=truefirings(2,:); truelabels=truefirings(3,:);
o.max_matching_offset = 2;    % in samples
[~,Q,acc,t] = times_labels_accuracy(truetimes,1+0*truelabels,times,1+0*times,o);
