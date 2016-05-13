% Compare sorters on Buzsaki.  Model for eventual compare_two_sorters
% Barnett 5/13/16

clear; dA = grab_buzsaki_dataset(2);

opts.detect_threshold = 4.0;    % stddev units
opts.freq_min = 400; opts.freq_max = 10000;  % don't make min too low
opts.samplerate = dA.samplerate;
opts.clip_size = 50;         % 2.5 ms
opts.detect_polarity = 'm';     % (simplesorter) prevents t-shifted duplicates
opts.sign=-1;
opts.detectability_threshold=3;  % why so low otherwise no events?? for (1)
opts.num_fea = 20;  % note only get 4-5 with nfea=10, but now get 5-6.
opts.artifacts = 1;
%opts.verb = 1;           % makes figure of clustering in PCA space (unordered)

% amazing how much variation there is in these sortings!

outdir = dA.outdir;
dA.outdir = [outdir 'A']; mkdir(dA.outdir);
dB = dA;
dB.outdir = [outdir 'B']; mkdir(dB.outdir);  % note distinct output directories
[dA.firings,infoA] = simplesorter(dA.timeseries,dA.outdir,opts);
[dB.firings,infoB] = jfm_april_sort(dB.timeseries,dB.outdir,opts);

% useful to standardize the ordering somewhat (rewrites firingsfile)...
[clipsA,templatesA] = reorderbytemplatenorm(dA.firings,infoA.prefile,opts);

dA.timeseries = [dA.outdir,'/pre.mda'];  % show pre-proc not raw
dB.timeseries = [dB.outdir,'/pre2.mda'];
o.verb = 1;
[fk Q perm iperm] = compare_two_sortings(dA,dB,o)

% reorder B's firingsfile...
firings = readmda(dB.firings); times=firings(2,:); labels=firings(3,:);
writemda64([0*times;times;perm(labels)],dB.firings);

% show sortings...
mv.raw=dA.timeseries; mv.pre=infoA.prefile; mv.filt=infoA.filtfile;
mv.firings=dA.firings; mv.samplerate=dA.samplerate; mountainview(mv);

mv.raw=dB.timeseries; mv.pre=infoB.prefile; mv.filt=infoB.filtfile;
mv.firings=dB.firings; mv.samplerate=dB.samplerate; mountainview(mv);
