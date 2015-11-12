function example1_mountainsort

%close all;

%Set some options and specify the input/output file names
basepath=[fileparts(mfilename('fullpath')),'/..'];
addpath(basepath);
opts.raw_dat=[basepath,'/../ms11d45A/ms11d45.dat'];
opts.channels=[37:52,68,69];
opts.adjacency_radius=2;
opts.raw_mda=[basepath,'/../ms11d45A/ms11d45A_pre.mda'];
opts.locations=[basepath,'/../ms11d45A/locations.mda'];
opts.num_cluster_features=6;
opts.clip_size=80;
opts.detection_threshold=5;
opts.detection_interval=40;
opts.detect_freq_min=600;
opts.detect_freq_max=2000;
opts.working_path=[basepath,'/../ms11d45A/working'];
opts.output_path=[basepath,'/../ms11d45A/output'];
opts.prewhiten=1;

%opts.timepoints=1:19e6; %Something seems to change around timepoint 19-20 million
%opts.timepoints=[1:2.9e6,3.0e6:19e6]; %Something seems to change around timepoint 19-20 million
%opts.timepoints=[1:2.9e6,3.0e6:5e6];
opts.timepoints=1:1e6;
%opts.timepoints=1:4.5e6;
%opts.timepoints=1:10e6;

if (~exist(opts.raw_dat,'file'))&&(~exist(opts.raw_mda,'file'))
    
    basepath_parent=fileparts(basepath);
    
    url='http://voms.simonsfoundation.org:50013/rXcGof5b0pDR4zkEDCevBLfo95RB7/frank_ucsf_example1/ms11d45A_pre.mda';
    
    fprintf('\n');
    fprintf('It appears that you do not have the raw data file on your computer.\n');
    
    fprintf('You can download the file from here: %s\n',url);
    fprintf('Then copy the file to: %s/ms11d45A/ms11d45A_pre.mda\n\n',basepath_parent);
    
    fprintf('Easier: You could do the following if you had curl:\n\n');
    fprintf('> mkdir %s/ms11d45A\n',basepath_parent);
    fprintf('> curl %s > %s/ms11d45A/ms11d45A_pre.mda\n\n',url,basepath_parent);
    
    fprintf('Once you have obtained the file (around 4GB), re-run this script.\n\n');
    
    warning('Raw file not found. See the above message.');
    return;
end;

mountainsort(opts);

end



