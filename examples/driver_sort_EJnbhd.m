% Show how to call the simple sorter on real EJ data (7-channel neighborhood)
% Barnett 3/21/16 for ISOSPLIT paper w/ Magland
% I guess can swap in different clusterings here? (or use those in validspike)
%
% If you want ground-truthed data from the same waveforms, see
%  driver_simplesorter.m

clear
d = grab_EJnbhd_dataset;   % get struct pointing to tmp files

opts.detect_threshold = 120;    % abs units (uV), not sure why
opts.freq_min = 100; opts.freq_max = inf;  % I think EJ already filtered anyway
opts.samplerate = d.samplerate;
opts.clip_size = 50;         % 2.5 ms
%opts.detect_meth = 'x';  % 'x' = don't do PCA-alignment
opts.detect_polarity = 'm';     % needed for EJ to prevent t-shifted duplicates
%opts.verb = 1;           % makes figure of clustering in PCA space (unordered)

% go into this code to swap clustering algs...
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

% [there is a duplicated type which is weird - is detection correct?
%  FIXED: polarity negative only fixed it]
