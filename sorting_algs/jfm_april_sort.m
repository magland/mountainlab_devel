function [firingsfile,info] = jfm_april_sort(tsfile,path,o)
% JFM_APRIL_SORT   JFM 4/8/16 sorter w/ std MS MDA interface.
%
% [firingsfile,info] = jfm_april_sort(rawfile,output_dir,o)
%
% Inputs:
%    rawfile - path to .mda of MxN raw signal data
%    output_dir - path to existing DIRECTORY where all output will be written
%    o - (optional) sorting options, including:
%        clip_size
%        detect_threshold (in sigma)
%        samplerate
%        detect_interval
%        min_shell_size
%        shell_increment
%        sign = 0,-1,+1
%        consolidation_factor - between 0 and 1
%        whiten = 0 or 1
%        fit = 0 or 1
%        artifacts = 0 or 1
%        merge_threshold, eg 0.9
%        outlier_threshold, eg 5
%        detectability_threshold, eg 5
%
% Outputs:
%    firingsfile - path to the firings.mda output file
%    info - struct with fields:
%           filtfile - filtered timeseries
%           prefile - path to the preprocessed timeseries (filt and whitened)
%
% Other m-files required: mscmd_*
%
% wrapped from jfm, Barnett 4/8/16. todo: break out some more opts

defo.clip_size = 60;       % default opts
defo.detect_threshold = 3.0;
defo.samplerate = 30000;
defo.detect_interval = 10;
defo.min_shell_size = 150;
defo.shell_increment = 1.5;
defo.sign = 0;
defo.consolidation_factor=0.9; %added 4/12/16
defo.whiten=1; %added 4/13/16
defo.fit=1; %added 4/13/16
defo.artifacts=0; %added 4/13/16
defo.merge_threshold=0.9; %added 4/13/16
defo.detectability_threshold=6; %added 4/13/16
defo.outlier_threshold=5; %added 4/13/16
if nargin<3 o=struct; end; o = ms_set_default_opts(o,defo); % setup opts

disp('sorter opts:'); o

o_filter.samplerate=o.samplerate;
o_filter.freq_min=300;
o_filter.freq_max=10000;

o_mask_out_artifacts.threshold=3;
o_mask_out_artifacts.interval_size=200;

o_detect.detect_threshold=o.detect_threshold;
o_detect.detect_interval=o.detect_interval;
o_detect.clip_size=o.clip_size;
o_detect.sign=o.sign;
o_detect.upsampling_factor = 5;   % jfm suggestion for detect3
o_detect.num_pca_denoise_components = 5;   % "
o_detect.pca_denoise_jiggle = 2;           % "

o_branch_cluster.clip_size=o.clip_size;
o_branch_cluster.min_shell_size=o.min_shell_size;
o_branch_cluster.shell_increment=o.shell_increment;
o_branch_cluster.num_features=3;
o_branch_cluster.detect_interval=o.detect_interval;
o_branch_cluster.consolidation_factor=o.consolidation_factor; %added 4/12/16

o_compute_detectability_scores.clip_size=o.clip_size;
o_compute_detectability_scores.shell_increment=o.shell_increment;
o_compute_detectability_scores.min_shell_size=o.min_shell_size;

o_compute_outlier_scores.clip_size=o.clip_size;
o_compute_outlier_scores.shell_increment=o.shell_increment;
o_compute_outlier_scores.min_shell_size=o.min_shell_size;

o_fit_stage.clip_size=o.clip_size;
o_fit_stage.min_shell_size=o.min_shell_size;
o_fit_stage.shell_increment=o.shell_increment;

o_merge_labels.clip_size=o.clip_size;
o_merge_labels.merge_threshold=o.merge_threshold;

o_filter_events.detectability_threshold=o.detectability_threshold;
o_filter_events.outlier_threshold=o.outlier_threshold;

tA=tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_bandpass_filter(tsfile,[path,'/pre1.mda'],o_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if o.artifacts
    mscmd_mask_out_artifacts([path,'/pre1.mda'],[path,'/pre1b.mda'],o_mask_out_artifacts);
else
    mscmd_copy([path,'/pre1.mda'],[path,'/pre1b.mda']);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if o.whiten
    mscmd_whiten([path,'/pre1b.mda'],[path,'/pre2.mda']);
else
    mscmd_copy([path,'/pre1b.mda'],[path,'/pre2.mda']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_detect3([path,'/pre2.mda'],[path,'/detect.mda'],o_detect);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_branch_cluster_v2([path,'/pre2.mda'],[path,'/detect.mda'],'',[path,'/firings1.mda'],o_branch_cluster);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_copy([path,'/firings1.mda'],[path,'/firings2.mda']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if o.fit
    mscmd_fit_stage([path,'/pre2.mda'],[path,'/firings2.mda'],[path,'/firings3.mda'],o_fit_stage);
else
    mscmd_copy([path,'/firings2.mda'],[path,'/firings3.mda']);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_merge_labels([path,'/pre2.mda'],[path,'/firings3.mda'],[path,'/firings4.mda'],o_merge_labels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_compute_outlier_scores([path,'/pre2.mda'],[path,'/firings4.mda'],[path,'/firings5.mda'],o_compute_outlier_scores);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_compute_detectability_scores([path,'/pre2.mda'],[path,'/firings5.mda'],[path,'/firings6.mda'],o_compute_detectability_scores);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_filter_events([path,'/firings6.mda'],[path,'/firings7.mda'],o_filter_events);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mscmd_copy([path,'/firings7.mda'],[path,'/firings.mda']);

info.filtfile = [path,'/pre1b.mda'];
info.prefile = [path,'/pre2.mda'];
firingsfile = [path,'/firings.mda'];

%mscmd_mda2txt([path,'/firings.mda'],[path,'/firings.txt']);

fprintf('Total time for sorting: %g sec\n',toc(tA));

%mv.mode='overview2';
%mv.raw=[path,'/pre0.mda'];
%mv.filt=[path,'/pre1b.mda'];
%mv.pre=[path,'/pre2.mda'];
%mv.firings=[path,'/firings.mda'];
%mv.samplerate=samplerate;
%mountainview(mv);
