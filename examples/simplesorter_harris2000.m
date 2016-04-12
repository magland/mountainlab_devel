% Show how to call the simple sorter on Harris 2000 tetrode.
% Barnett 4/11/16

clear; d = grab_harris2000_dataset;

opts.samplerate = d.samplerate;
opts.clip_size = 50;
opts.detect_polarity = 'm';   % needed for good detection
[firingsfile,info] = simplesorter(d.timeseries,d.outdir,opts);

% useful to standardize the ordering somewhat...
[clips,templates] = reorderbytemplatenorm(firingsfile,info.prefile,opts);

% load and view the EC input signal with firings output, and waveforms
d.timeseries = info.prefile;   % switch to preproc timeseries
Y = readmda(d.timeseries);
firings = readmda(firingsfile); times=firings(2,:); labels=firings(3,:);
figure; ms_view_clusters(ms_event_features(clips,3,opts),labels); % prepr-clips 
spikespy({Y,times,labels,'simple sorter'});
K = max(labels); pops = histc(labels,1:K); disp('populations n_k vs k:');
fprintf('\t%d',1:K); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');
clips = ms_extract_clips2(Y,times,opts.clip_size);  % raw clips & wf's from them
W = ms_templates(clips,labels);
figure; ms_view_templates(W,struct('showcenter',1)); title('sorted W from raw');

mv.mode='overview2'; mv.raw=d.timeseries; mv.pre=info.prefile;  % nice view...
mv.firings=firingsfile; mv.samplerate=d.samplerate; mountainview(mv);

%o.dtau = 1e-3*opts.samplerate; o.taumax = 20e-3*opts.samplerate;
%show_crosscorr(labels,times,[],o);  % seems way too slow

% set up 2 datasets to do comparison to IC truth...
da = d; da.firings = d.truefirings; da.name = [d.name ' IC'];
db = d; db.firings = firingsfile; db.name = [d.name ' sorted'];
[fk Q perm] = compare_two_sortings(da,db);
fprintf('Accuracies vs label k are... \n')
fprintf('k   :'); fprintf('\t%d',1:numel(fk)), fprintf('\n');
fprintf('f_k :'); fprintf('\t%.3f',fk), fprintf('\n');
