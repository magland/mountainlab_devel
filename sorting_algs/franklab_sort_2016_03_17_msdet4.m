function [firings,info] = franklab_sort_2016_03_17_msdet4(signal,outdir,o)
% matlab driver for jfm's franklab sorter from March 17 2016, with detection
%  updated to mscmd_detect3 - this is in progress.
%
% [firingsfile,info] = franklab_sort_2016_03_17(rawfile,output_dir,opts)
%
% Inputs:
%    rawfile - path to .mda of MxN raw signal data
%    output_dir - path to existing DIRECTORY where all output will be written
%    sort_opts - (optional) sorting options, see def_sort_opts in this
%                   script
%
% Outputs:
%    firingsfile - path to the firings.mda output file
%    info - struct with fields:
%           prefile - path to the preprocessed raw data file
%
% todo: break out some hardwired opts

mfile_path=fileparts(mfilename('fullpath'));
path=outdir;

clip_size=100;
detect_interval=10;
detectability_threshold=4;
shell_increment=1.5;
min_shell_size=150;
samplerate=o.samplerate;
o_filter.samplerate=samplerate;
o_filter.freq_min=300;
o_filter.freq_max=6000;
o_mask_out_artifacts.threshold=3;
o_mask_out_artifacts.interval_size=100;

o_detect.detect_threshold=3.5;       % opts now for ms_detect4
o_detect.detect_interval=detect_interval;
o_detect.clip_size=clip_size;
o_detect.polarity = 'm';
o_detect.beta = 10;
o_detect.num_features = 5;    % since individ chan
o_detect.individual_channels = 1;        % for jfm

o_branch_cluster.clip_size=clip_size;
o_branch_cluster.min_shell_size=min_shell_size;
o_branch_cluster.shell_increment=shell_increment;
o_branch_cluster.num_features=10;
o_branch_cluster.detect_interval=detect_interval;
o_remove_duplicate_clusters.max_dt=6;
o_remove_duplicate_clusters.overlap_threshold=0.25;
o_remove_noise_subclusters.clip_size=clip_size;
o_remove_noise_subclusters.detectability_threshold=detectability_threshold;
o_remove_noise_subclusters.shell_increment=shell_increment;
o_remove_noise_subclusters.min_shell_size=min_shell_size;
o_compute_outlier_scores.clip_size=clip_size;
o_compute_outlier_scores.shell_increment=shell_increment;
o_compute_outlier_scores.min_shell_size=min_shell_size;

tA=tic;

mscmd_bandpass_filter(signal,[path,'/pre1.mda'],o_filter);
mscmd_mask_out_artifacts([path,'/pre1.mda'],[path,'/pre1b.mda'],o_mask_out_artifacts);
mscmd_whiten([path,'/pre1b.mda'],[path,'/pre2.mda']);

disp('--- MS_DETECT4 ---'); % insert pure-Matlab detection...
[times,peakchans] = ms_detect4(readmda([path,'/pre2.mda']),o_detect);
writemda([peakchans;times], [path,'/detect.mda'], 'float64');  % just 2 rows

mscmd_branch_cluster_v2([path,'/pre2.mda'],[path,'/detect.mda'],'',[path,'/firings1.mda'],o_branch_cluster);
mscmd_remove_duplicate_clusters([path,'/pre2.mda'],[path,'/firings1.mda'],[path,'/firings2.mda'],o_remove_duplicate_clusters);  % ahb added 1st arg
mscmd_remove_noise_subclusters([path,'/pre2.mda'],[path,'/firings2.mda'],[path,'/firings3.mda'],o_remove_noise_subclusters);
mscmd_compute_outlier_scores([path,'/pre2.mda'],[path,'/firings3.mda'],[path,'/firings4.mda'],o_compute_outlier_scores);

mscmd_copy([path,'/firings4.mda'],[path,'/firings.mda']);
mscmd_mda2txt([path,'/firings.mda'],[path,'/firings.txt']);

fprintf('Total time for sorting: %g sec\n',toc(tA));

% output fnames...
firings = [path,'/firings.mda'];
info.prefile = [path,'/pre2.mda'];
info.filtfile = [path,'/pre1.mda'];

if 0
  mv.mode='overview2';
  mv.raw=[path,'/pre0.mda'];
  mv.filt=[path,'/pre1b.mda'];
  mv.pre=[path,'/pre2.mda'];
  mv.firings=[path,'/firings.mda'];
  mv.samplerate=samplerate;
  mountainview(mv);
end
