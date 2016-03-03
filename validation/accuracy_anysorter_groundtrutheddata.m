function fk = accuracy_anysorter_groundtrutheddata(sortfunc,datasetname,outdir,o_acc,o_sorter)
% ACCURACY_ANYSORTER_GROUNDTRUTHEDDATA  test accuracy of a sorter on synth data.
%
% fk = accuracy_anysorter_groundtrutheddata(sortfunc,datasetname,outdir,o_acc,
%      o_sorter)
% measures any sorter's accuracy vs known firing times on synthetic demo data.
% Also pops up a bunch of windows showing accuracy.
%
% Inputs:
%   sortfunc - function handle to a sorter with the interface
%                [firingsfile,prefile] = sorter(rawfile,outdir);
%
%  datasetname - text string passed to get_default_dataset
%  outputdir - working directory where MDA files will go
%  o_acc - accuracy measuring options:
%               o_acc.samplefreq - sample rate in Hz
%               o_acc.T          - clip size
%  o_sorter (optional) - passed to last arg of sorter, as in:
%                [firingsfile,prefile] = sorter(rawfile,outdir,opts);
% Outputs:
%  fk - (1xK) accuracy metrics on the neuron types labeled by ground truth labels
%
% Run without options runs a demo test

% Barnett 3/3/16 based on accuracy_sort_demotimeseries.m
if nargin==0, test_accuracy_anysorter_groundtrutheddata; return; end
if nargin<2 | isempty(datasetname), datasetname = 'EJ K7'; end
if nargin<3, o_acc=[]; end
if ~isfield(o_acc,'samplefreq'), o_acc.samplefreq = 2e4; end
if ~isfield(o_acc,'T'), o_acc.T = 50; end
if ~isfield(o_acc,'max_matching_offset'), o_acc.max_matching_offset = 10; end  % in samples; 0.5 ms
if nargin<4, o_sorter=[]; end

[Yfile truefiringsfile trueWfile] = get_default_dataset(datasetname);
if ~exist(outdir,'dir') mkdir(outdir); end

% load & view the data and truth about it...
Y=readmda(Yfile);                 % raw
truefirings=readmda(truefiringsfile);       % read sort output files
truetimes=truefirings(2,:); truelabels=truefirings(3,:);
clips=ms_extract_clips(Y,round(truetimes),o_acc.T);  % nearest integer t-pt
trueX = ms_event_features(clips,3);       % 3 features for vis
figure; set(gcf,'position',[1000 1000 1000 1000]); tsubplot(2,2,1);
ms_view_clusters(trueX,truelabels); title('true labels in fea space');
Wtrue = readmda(trueWfile);
tsubplot(2,2,2); ms_view_templates(Wtrue), title('true W'); drawnow

firingsfile=[outdir '/firings.mda'];    % make MDA for sorter to read

[firingsfile,prefile] = sortfunc(Yfile,outdir,o_sorter);  % DO IT

firings=readmda(firingsfile);       % read sort output file
times=firings(2,:); labels=firings(3,:);
clips=ms_extract_clips(Y,times,o_acc.T);
X = ms_event_features(clips,3);       % 3 features, just for viz

if 0 %%%% View output
  mv.mode='overview2'; mv.raw=prefile; mv.firings=firingsfile;
  mv.sampling_freq = o_acc.samplefreq;
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
if 1     % old matlab validspike way
  [perm Q acc t] = times_labels_accuracy(truetimes,truelabels,times,labels,o_acc);
else     % C fast way, but not full times-labels best search
  mscmd_confusion_matrix(truefiringsfile,firingsfile,[outdir '/accconfmat.mda'],o_acc.max_matching_offset);
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

%addpath ~/spikespy/matlab/  % prefer old spikespy
spikespy({Y,truetimes,truelabels,'Y, truth'},{Y,times,perm(labels),'Y, sorted'},{[t.tmiss';t.tfals';t.twrng'],[1+0*t.tmiss';2+0*t.tfals';3+0*t.twrng']});
%rmpath ~/spikespy/matlab/

%%%%%%%%%%%%%%%%%%%%%%%%


function test_accuracy_anysorter_groundtrutheddata

mfile_path=fileparts(mfilename('fullpath'));
addpath([mfile_path,'/../demo/demo_sort_001']);
fk = accuracy_anysorter_groundtrutheddata(@ds001_sort,[],'/tmp/output',[],[]);
fprintf('done; fk accuracies are... \n')
fk
