% Show how to call the simple sorter on Buzsaki.
% Barnett 3/24/16

clear
d = grab_buzsaki_dataset(2);

opts.detect_threshold = 3.5;    % stddev units
opts.freq_min = 400; opts.freq_max = 10000;  % don't make min too low
opts.samplerate = d.samplerate;
opts.clip_size = 50;         % 2.5 ms
opts.detect_polarity = 'm';     % (simplesorter) prevents t-shifted duplicates
opts.sign=-1;
opts.detectability_threshold=0;  % why so low otherwise no events?? for (1)
%opts.verb = 1;           % makes figure of clustering in PCA space (unordered)

% go into this code to swap clustering algs...
[firingsfile,info] = simplesorter(d.timeseries,d.outdir,opts);
%[firingsfile,info] = jfm_april_sort(d.timeseries,d.outdir,opts);

% useful to standardize the ordering somewhat (rewrites firingsfile)...
[clips,templates] = reorderbytemplatenorm(firingsfile,info.prefile,opts);

% load and view the filtered signal with firings output, and waveforms
Y = readmda(info.prefile);
firings = readmda(firingsfile); times=firings(2,:); labels=firings(3,:);
figure; ms_view_clusters(ms_event_features(clips,3,opts),labels); % prepr-clips 
spikespy({Y,times,labels,'simple sorter'});
K = max(labels); pops = histc(labels,1:K); disp('populations n_k vs k:');
fprintf('\t%d',1:K); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');
%clips = ms_extract_clips2(Y,times,opts.clip_size);  % raw clips & wf's from them
%W = ms_templates(clips,labels);
figure; ms_view_templates(templates,struct('showcenter',1,'pops',pops));title('sorted W from pre');

mv.raw=d.timeseries; mv.pre=info.prefile; mv.filt=info.filtfile; % nice view...
mv.firings=firingsfile; mv.samplerate=d.samplerate; mountainview(mv);
