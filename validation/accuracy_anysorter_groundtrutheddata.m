function fk = accuracy_anysorter_groundtrutheddata(sortfunc,dataset,o_acc,o_sorter)
% ACCURACY_ANYSORTER_GROUNDTRUTHEDDATA  meas accuracy of sorter on gnd-truth data
%
% fk = accuracy_anysorter_groundtrutheddata(sortfunc,dataset,o_acc,o_sorter)
%   measures any sorter's accuracy vs known firing times on synthetic demo data.
%   Also pops up a bunch of windows showing accuracy.
%
% Inputs:
%  sortfunc - function handle to a sorter with the interface
%                [firingsfile,prefile] = sorter(rawfile,outdir);
%  dataset - struct with fields giving dataset with its ground truth:
%            dataset.signal - MDA filename or M*N array, raw EC signal data
%            dataset.truefirings - MDA filename or 4*Ns array, true firings
%            dataset.truewaveforms - MDA filename or M*T*K array, true waveforms
%                                    (optional)
%            dataset.samplerate - samplerate in samples/sec for raw data
%            dataset.outdir - path to existing directory for MDA output files
%  o_acc - accuracy measuring options:
%               o_acc.T          - clip size
%  o_sorter (optional) - passed to last argument of sortfunc, as in:
%                [firingsfile,prefile] = sortfunc(rawfile,outdir,o_sorter)
%                Note: the sorter opts will already be passed the samplerate
%                      as given in the dataset struct.
% Outputs:
%  fk - (1xK) accuracy metrics on the neuron types labeled by ground truth labels
%
% Run without options runs a test on demo data

% Barnett 3/3/16 based on accuracy_sort_demotimeseries.m
% 3/16/16 Barnett changed to use non-integer firing time resampling
% todo: * make various switches below into options.
% 3/18/16 struct not text input for dataset. truewaveforms optional

if nargin==0, accuracy_simplesorter; return; end  % is in validation/
if nargin<2 | isempty(dataset), dataset = demo_dataset; end
if nargin<3, o_acc=[]; end
if ~isfield(o_acc,'T'), o_acc.T = 50; end
if ~isfield(o_acc,'max_matching_offset'), o_acc.max_matching_offset = 10; end  % in samples; eg 0.5 ms at 20 kHz
if nargin<4, o_sorter=[]; end

outdir = dataset.outdir;
samplerate = dataset.samplerate;
o_sorter.samplerate = samplerate;
o_sorter.samplefreq = samplerate;   % remove when that fieldname obsolete

% load & view the data and truth about it...
Y = arrayify(dataset.signal);                    % raw
truefirings=arrayify(dataset.truefirings);       % read sort output files
truetimes=truefirings(2,:); truelabels=truefirings(3,:);
betaonesnap = 10;           % upsampling for extract clips
clips=ms_extract_clips2(Y,truetimes,o_acc.T,[],betaonesnap);  % real t, resamples
%clips=ms_extract_clips(Y,round(truetimes),o_acc.T);  % int t
trueX = ms_event_features(clips,3);       % 3 features for vis
figure; set(gcf,'position',[1000 1000 1000 1000]); tsubplot(2,2,1);
ms_view_clusters(trueX,truelabels); title('true labels in fea space');
Wopts.showcenter = 1;
if isfield(dataset,'truewaveforms')
  Wtrue = arrayify(dataset.truewaveforms);
  tsubplot(2,2,2); ms_view_templates(Wtrue,Wopts)
  title('true W'); drawnow
end
trueK = max(truelabels); pops = histc(truelabels,1:trueK);
fprintf('true populations for each label:\n');
fprintf('\t%d',1:trueK); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');

firingsfile=[outdir '/firings.mda'];    % make MDA for sorter to read

disp('call sorter:')
[firingsfile,prefile] = sortfunc(fnameify32(dataset.signal,outdir),outdir,o_sorter);  % DO IT

firings=readmda(firingsfile);       % read sort output file
times=firings(2,:); labels=firings(3,:);
clips=ms_extract_clips2(Y,times,o_acc.T,[],betaonesnap);
X = ms_event_features(clips,3);       % 3 features, just for viz

if 0 %%%% View output
  mv.mode='overview2'; mv.raw=prefile; mv.firings=firingsfile;
  mv.sampling_freq = samplefreq;
  ms_mountainview(mv);
end

templates = ms_templates(clips,labels);    % get mean waveforms

K = max(labels);  % # claimed neurons
ii = labels>0;  % decide handling unclassified
if 0                % kill unclassified
  times = times(ii); labels=labels(ii);
else                % treat unclass as a new type
  labels(~ii) = K+1; K = max(labels);
end

fprintf('compute accuracy vs known ground-truth...\n')
if 1     % old matlab validspike way (3-pass)
  [perm Q acc t]=times_labels_accuracy(truetimes,truelabels,times,labels,o_acc);
else     % C fast way, but doesn't do full times-labels best search
  mscmd_confusion_matrix(truefiringsfile,firingsfile,...
                         [outdir '/accconfmat.mda'],o_acc.max_matching_offset);
  Q = readmda([outdir '/accconfmat.mda']);   % un-permed
  [perm Q] = bestcolpermconfmat(Q);
  disp('best permed extended confusion matrix:'), Q
  [~,acc] = labels_similarity(Q);
end
fk = acc.p;   % our accuracy measure

pops = histc(perm(labels),1:K);
fprintf('sorter populations for each label (best permuted):\n');
fprintf('\t%d',1:K); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');

tsubplot(2,2,3); ms_view_clusters(X,perm(labels));  % show accuracy
title('sorted best-permed labels');
[~,iperm] = sort(perm(1:K));
tsubplot(2,2,4); ms_view_templates(templates(:,:,iperm),Wopts);
title('best-permed W from sorting');

figure; set(gcf,'position',[1000 500 1000 400]);
subplot(1,2,1); imagesc(Q); colorbar;ylabel('true label');xlabel('sorted label');
hold on; plot([.5,K+.5;K+1.5,K+.5], [trueK+.5,.5;trueK+.5,K+1.5],'w-');
title('best extended accuracy confusion matrix');
subplot(1,2,2); plot(fk,'.','markersize',20); axis([1 K 0 1]);
xlabel('k (best permuted output label)'); ylabel('accuracy metric f_k');

%addpath ~/spikespy/matlab/  % prefer old spikespy
spikespy({Y,truetimes,truelabels,'Y, truth'},{Y,times,perm(labels),'Y, sorted'},{[t.tmiss';t.tfals';t.twrng'],[1+0*t.tmiss';2+0*t.tfals';3+0*t.twrng']});
%rmpath ~/spikespy/matlab/

%%%%%%%%%%%%%%%%%%%%%%%%


% helpers...

function fname = fnameify32(X,outdir)
% FNAMEIFY  if array, writes to file and returns filename, otherwise keeps name

% v crude for now.

if ischar(X) || isstring(X)
  fname = X;
else
  fname = [outdir,'/',num2str(randi(1e10)),'.mda'];  % random filename
  writemda32(X,fname);
end
%%%%%%%%

function fname = fnameify64(X,outdir)
% FNAMEIFY  if array, writes to file and returns filename, otherwise keeps name

% v crude for now.

if ischar(X) || isstring(X)
  fname = X;
else
  fname = [outdir,'/',num2str(randi(1e10)),'.mda'];  % random filename
  writemda64(X,fname);
end
%%%%%%%%

function X = arrayify(X)
if ischar(X) || isstring(X), X = readmda(X); end
%%%%%
