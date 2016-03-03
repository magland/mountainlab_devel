function accuracy_sort_demotimeseries
% measure accuracy vs known firing times on fixed synthetic demo data.
% Currently uses some C executables in the chain, shuffling data to MDA & back.
% Barnett 2/25/16. 2/26/16 playing w/ sorting chain & viz.
% 3/3/16 J's ds001

%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'));
addpath([mfile_path,'/../franklab/sort_002_multichannel']);
addpath([mfile_path,'/../demo/demo_sort_001']);
rng(1);  % fix the seed
[Yfile truefiringsfile trueWfile o.samplefreq] = get_default_dataset('EJ K7');
outdir=[mfile_path,'/output'];
if ~exist(outdir,'dir') mkdir(outdir); end;

% load & view the data and truth about it...
Y=readmda(Yfile);                 % raw
truefirings=readmda(truefiringsfile);       % read sort output files
truetimes=truefirings(2,:); truelabels=truefirings(3,:);
T = 50;          % clip size, 2.5 ms (truth at most 46)
clips=ms_extract_clips(Y,round(truetimes),T);  % nearest integer t-pt
trueX = ms_event_features(clips,3);       % 3 features for vis
figure; set(gcf,'position',[1000 1000 1000 1000]); tsubplot(2,2,1);
ms_view_clusters(trueX,truelabels); title('true labels in fea space');
Wtrue = readmda(trueWfile);
tsubplot(2,2,2); ms_view_templates(Wtrue), title('true W'); drawnow

% common sorting opts
o.freq_min = 100; o.freq_max = 5000;   % filter for det, not clustering

if 1         % Standard simple matlab chain
  o.clip_size=T; o.detect_interval = 20;
  o.detect_threshold = 80;     % absolute (uV);  80 uV collects a noise clus
  %o.detect_threshold = 3.5;   % only if whitened, in sigma units
  %Y = ms_filter(Y,o); Y = ms_whiten(Y);  % if want match what ds001 does
  Yf = ms_filter(Y,o);
  if 0   % the C interface to detect
    prefile = [outdir,'/pre.mda'];
    writemda(Yf,prefile);          % in case Y differs from raw
    mscmd_detect(prefile,[outdir,'/detect.mda'],o);
    times=readmda([outdir,'/detect.mda']);
    times=times(2,:);
  else
    times = ms_detect(Yf,o);
  end
  if 0      % the C interface to get clips
    writemda([0*times; times], [outdir,'/detect.mda']);
    clip_opts.clip_size=T;
    mscmd_extract_clips(Yfile,[outdir,'/detect.mda'],[outdir,'/clips.mda'],clip_opts);
    clips=readmda([outdir,'/clips.mda']);
  else      % easier interface (pure Matlab)
    clips = ms_extract_clips(Y,times,T);
  end
  numfea = 12; X = ms_event_features(clips,numfea);
  [labels,info] = isosplit2(X);           % cluster
  firingsfile=[outdir '/firings.mda'];    % make MDA for MV and mscmd_* to read
  writemda(TL2F(times,labels),firingsfile);
  prefile = Yfile;    % since no preproc
  
elseif 1     % J's ds001 chain
  o.detect_threshold=3;
  o.detect_interval=20;   % in samples
  o.shell_increment=1.0;    % doubled so faster
  firingsfile=[outdir '/firings.mda'];    % make MDA for MV and mscmd_* to read
  [firingsfile,prefile]=ds001_sort(Yfile,outdir,o);
  firings=readmda(firingsfile);       % read sort output file
  times=firings(2,:); labels=firings(3,:);
  clips=ms_extract_clips(Y,times,T);
  X = ms_event_features(clips,3);       % 3 features, just for vis
  
else          % J's experimental stuff
  sort_opts=struct;
  sort_opts.detect.detect_threshold=4;   % in sigmas; don't vary much
  sort_opts.detect.detect_interval=15;
  sort_opts.detectability_threshold=0; % controls exclusion of noise clusters
  sort_opts.plausibility_threshold=inf;  % exc
  sort_opts.clip_size=T;
  [firingsfile,prefile]=sort_002_multichannel(Yfile,outdir,sort_opts);
  firings=readmda(firingsfile);       % read sort output files
  times=firings(2,:); labels=firings(3,:);
  clips=ms_extract_clips(Y,times,T);
  X = ms_event_features(clips,3);       % 3 features, just for vis
end

if 0 %%%% View output
  mv.mode='overview2';
  mv.raw=prefile;
  mv.firings=firingsfile;
  ms_mountainview(mv);
end

templates=ms_templates(clips,labels);    % mean waveforms

K = max(labels);  % # claimed neurons
ii = labels>0;  % decide handling unclassified
if 0                % kill unclassified
  times = times(ii); labels=labels(ii);
else                % treat unclass as a new type
  labels(~ii) = K+1; K = max(labels);
end

fprintf('compute accuracy vs known ground-truth...\n')
o.max_matching_offset = 10;    % in samples; 0.5 ms
if 1     % old matlab validspike way
  [perm Q acc t] = times_labels_accuracy(truetimes,truelabels,times,labels,o);
else     % C fast way, but not full times-labels best search
  mscmd_confusion_matrix(truefiringsfile,firingsfile,[outdir '/accconfmat.mda'],o.max_matching_offset);
  Q = readmda([outdir '/accconfmat.mda']);   % un-permed
  [perm Q] = bestcolpermconfmat(Q);
  disp('best permed extended confusion matrix:'), Q
  [~,acc] = labels_similarity(Q);
end
fk = acc.p;   % our accuracy measure

tsubplot(2,2,3); ms_view_clusters(X,perm(labels));  % show accuracy
title('sorted best-permed labels');
[~,iperm] = sort(perm(1:K));
tsubplot(2,2,4); ms_view_templates(templates(:,:,iperm));
title('best-permed W from sorting');

figure; set(gcf,'position',[1000 500 1000 400]);
subplot(1,2,1); imagesc(Q); colorbar; ylabel('truth'); xlabel('sorted');
title('best extended accuracy confusion matrix');
subplot(1,2,2); plot(fk,'.','markersize',20); axis([1 K 0 1]);
xlabel('k'); ylabel('acc f_k');

addpath ~/spikespy/matlab/  % prefer old spikespy
spikespy({Y,truetimes,truelabels,'Y, truth'},{Y,times,perm(labels),'Y, sorted'},{[t.tmiss';t.tfals';t.twrng'],[1+0*t.tmiss';2+0*t.tfals';3+0*t.twrng']});
rmpath ~/spikespy/matlab/

%keyboard    % use dbcont to continue
end
