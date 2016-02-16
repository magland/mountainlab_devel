% getting to know MS via spikesort EJ data. Run from MS top directory.
% Barnett 2/16/16
clear; verb = 1;  % verbosity (0=no plots, etc)x

% convert EJ data format to MDA:
rawpath='../datasets/EJ';
path0 = '../datasets/EJ/MSdata';     % MS output sandbox for datasets
if ~exist(path0,'dir'), system(['mkdir ',path0]); end
head = '2005-04-26_elec359';          % filehead
if ~exist([path0,'/pre0.mda'],'file')
  info = extract_raw_EJ(head,rawpath,path0);
end
  
o_filter.samplefreq = readmda([path0,'/samplerate.mda']);
o_filter.freq_min=300;
o_filter.freq_max=inf;
mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre1.mda'],o_filter);
%if verb, Y = readmda([path0,'/pre1.mda']); spikespy(Y); end

% detect events as times on each channel...
o_detect.threshold = 100;        % absolute
o_detect.individual_channels = 0;  % max across all channels
mscmd_detect([path0,'/pre1.mda'],[path0,'/detect.mda'],o_detect);
CT = readmda([path0,'/detect.mda']); T=CT(2,:); clear CT
fprintf('detect found %d events\n',numel(T))
if verb
  Y = readmda([path0,'/pre1.mda']); spikespy({Y,T,0*T});  % label 0 is unclass
end

o_extract_clips.clip_size=60;
mscmd_extract_clips([path0,'/pre1.mda'],[path0,'/detect.mda'],[path0,'/clips.mda'],o_extract_clips);
X = readmda([path0,'/clips.mda']); spikespy(X);


