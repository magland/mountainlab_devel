function mountainsort_example1

%close all;

%Set some options and specify the input/output file names
basepath=fileparts(mfilename('fullpath'));
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

%opts.timepoints=1:19e6; %Something seems to change around timepoint 19-20 million
%opts.timepoints=[1:2.9e6,3.0e6:19e6]; %Something seems to change around timepoint 19-20 million
opts.timepoints=[1:2.9e6,3.0e6:5e6];
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

function step5_residual(opts,data)

timerA=tic;

fprintf('Step 5: Residual...\n');

consolidate_times_path=opts.consolidate_times_path;
consolidate_labels_path=opts.consolidate_labels_path;
consolidate_templates_path=opts.consolidate_templates_path;
consolidate_load_channels_path=opts.consolidate_load_channels_path;
residual_data_path=opts.residual_data_path;
residual_times_path=opts.residual_times_path;
residual_labels_path=opts.residual_labels_path;

AM=readmda(opts.adjacency);
M=size(AM,1);

times=readmda(consolidate_times_path);
labels=readmda(consolidate_labels_path);
templates=readmda(consolidate_templates_path);
load_channels=readmda(consolidate_load_channels_path);

fprintf('Forming residual...\n');
T=size(templates,2);
tt1=-ceil((T)/2);
tt2=tt1+T-1;
tt=tt1:tt2;
X=data.X;
for j=1:length(times)
    X(:,times(j)+tt)=X(:,times(j)+tt)-templates(:,:,labels(j));
end;

fprintf('Writing %s...\n',residual_data_path);
writemda(X,residual_data_path);

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');


end


