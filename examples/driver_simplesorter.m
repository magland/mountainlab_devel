% Show how to call the simple sorter on the demo data.
% Barnett 3/18/16

clear
d = demo_dataset;   % get struct pointing to demo data files

opts.samplerate = d.samplerate;
opts.clip_size = 50;
%opts.detect_polarity = 'm';
[firingsfile,info] = simplesorter(d.signal,d.outdir,opts);

% useful to standardize the ordering somewhat...
[clips,templates] = reorderbytemplatenorm(firingsfile,info.prefile,opts);

% load and view the EC input signal with firings output, and waveforms
Y = readmda(d.signal);
firings = readmda(firingsfile); times=firings(2,:); labels=firings(3,:);
figure; ms_view_clusters(ms_event_features(clips,3,opts),labels); % prepr-clips 
spikespy({Y,times,labels,'simple sorter'});
K = max(labels); pops = histc(labels,1:K); disp('populations n_k vs k:');
fprintf('\t%d',1:K); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');
clips = ms_extract_clips2(Y,times,opts.clip_size);  % raw clips & wf's from them
W = ms_templates(clips,labels);
figure; ms_view_templates(W,struct('showcenter',1)); title('sorted W from raw');

mv.mode='overview2'; mv.raw=d.signal; mv.pre=info.prefile;  % nice view...
mv.firings=firingsfile; mv.samplerate=d.samplerate; mountainview(mv);
