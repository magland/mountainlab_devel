function fk = accuracy_anysorter_groundtrutheddata(sortfunc,dataset,outdir,o_acc,o_sorter)
% ACCURACY_ANYSORTER_GROUNDTRUTHEDDATA  test accuracy of a sorter on synth data.
%
% fk = accuracy_anysorter_groundtrutheddata(sortfunc,dataset,outdir,o_acc,
%      o_sorter)
% measures any sorter's accuracy vs known firing times on synthetic demo data.
% Also pops up a bunch of windows showing accuracy.
%
% Inputs:
%   sortfunc - function handle to a sorter with the interface
%                [firingsfile,prefile] = sorter(rawfile,outdir);
%  dataset - struct with fields giving names of dataset with its ground truth:
%            dataset.raw - filename for raw data MDA
%            dataset.truefirings - filename for true firings MDA
%            dataset.truewaveforms - filename for true waveforms MDA
%            dataset.samplerate - samplerate in Hz for raw data
%            OR, if dataset is a text string, is passed to get_default_dataset
%            and chooses the name of generated or stored ground-truthed data.
%  outputdir - working directory where MDA files will go
%  o_acc - accuracy measuring options:
%               o_acc.T          - clip size
%  o_sorter (optional) - passed to last arg of sorter, as in:
%                [firingsfile,prefile] = sorter(rawfile,outdir,opts)
%                Note: the sorter opts will also be passed samplerate
%                      as given by the dataset.
% Outputs:
%  fk - (1xK) accuracy metrics on the neuron types labeled by ground truth labels
%
% Run without options runs a demo test

% Barnett 3/3/16 based on accuracy_sort_demotimeseries.m
% 3/16/16 Barnett changed to use non-integer firing time resampling
% todo: * make various switches below into options.
% 3/17/16 - struct input for dataset
if nargin==0, test_accuracy_anysorter_groundtrutheddata; return; end
if nargin<2 | isempty(dataset), dataset = 'EJ K7'; end
if nargin<3, o_acc=[]; end
if ~isfield(o_acc,'T'), o_acc.T = 50; end
if ~isfield(o_acc,'max_matching_offset'), o_acc.max_matching_offset = 10; end  % in samples; 0.5 ms
if nargin<4, o_sorter=[]; end

if isstruct(dataset)   % pull fields
  Yfile = dataset.raw;
  truefiringsfile = dataset.truefirings;
  trueWfile = dataset.truewaveforms;
  samplerate = dataset.samplerate;
else
  [Yfile truefiringsfile trueWfile samplerate] = get_default_dataset(dataset);
end
if ~isfield(o_sorter,'samplerate'), o_sorter.samplerate = samplerate;
  o_sorter.samplefreq = samplerate;   % remove when that fieldname obsolete
end
if ~exist(outdir,'dir') mkdir(outdir); end

% load & view the data and truth about it...
Y=readmda(Yfile);                 % raw
truefirings=readmda(truefiringsfile);       % read sort output files
truetimes=truefirings(2,:); truelabels=truefirings(3,:);
betaonesnap = 10;           % upsampling for extract clips
clips=ms_extract_clips2(Y,truetimes,o_acc.T,[],betaonesnap);  % real t, resamples
%clips=ms_extract_clips(Y,round(truetimes),o_acc.T);  % int t
trueX = ms_event_features(clips,3);       % 3 features for vis
figure; set(gcf,'position',[1000 1000 1000 1000]); tsubplot(2,2,1);
ms_view_clusters(trueX,truelabels); title('true labels in fea space');
Wtrue = readmda(trueWfile);
tsubplot(2,2,2); o.showcenter=1; ms_view_templates(Wtrue,o)
title('true W'); drawnow
trueK = max(truelabels); pops = histc(truelabels,1:trueK);
fprintf('true populations for each label:\n');
fprintf('\t%d',1:trueK); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');

firingsfile=[outdir '/firings.mda'];    % make MDA for sorter to read

disp('call sorter:')
[firingsfile,prefile] = sortfunc(Yfile,outdir,o_sorter);  % DO IT

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
tsubplot(2,2,4); ms_view_templates(templates(:,:,iperm),o);
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


function test_accuracy_anysorter_groundtrutheddata

mfile_path=fileparts(mfilename('fullpath'));
addpath([mfile_path,'/../demo/demo_sort_001']);
fk = accuracy_anysorter_groundtrutheddata(@ds001_sort,[],'/tmp/output',[],[]);
fprintf('done. fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');
