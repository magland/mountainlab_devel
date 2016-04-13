function [fk Q perm] = compare_two_sortings(da,db,o)
% COMPARE_TWO_SORTINGS  accuracy metrics between two sortings treating one true
%
% fk = compare_two_sortings(da,db,o) compares sorted output in db treating that
%  in da as truth. Produces several figures and text output. Dataset B may
%  have unclassified labels, which are treated as a new label type.
%
% Inputs:
%  da, db - each a dataset struct with at least the following fields:
%           timeseries - MDA filename or M*N array, raw EC signal
%           firings    - MDA filename or 4*Ns array, firings
%                        (row 1 is peak channels, row 2 firing times t_j, row 3
%                        is firing identities k_j, row 4 is firing amplitudes)
%           samplerate - timeseries samples/sec (at least db must have this).
%           name       - optional string giving short name of sortings
%  o - optional struct with optional fields:
%      T - clip size in samples
%      max_matching_offset - passed to accuracy measuring
%      betaonesnap - passed to clip extraction
%      verb - verbosity: 0 (text only), 1 (figure), 2 (figs+spikespy)
% Output:
%  fk - (1xK) accuracy metrics of d2 treating d1 as ground truth.
%  Q - extended best-permuted confusion matrix
%  perm - best permutation of B's labels to match A
%
% Without arguments, does self-test
%
% Used by: accuracy_anysorter_groundtrutheddata.m

% todo: update mscmd_conf_mat to consider better perms
%
% Barnett 4/5/16

if nargin==0, test_compare_two_sortings; return; end
if nargin<3, o = []; end
if ~isfield(o,'verb'), o.verb = 1; end
if ~isfield(o,'T'), o.T = 50; end
if ~isfield(o,'max_matching_offset'), o.max_matching_offset = 10; end  % in samples; eg 0.5 ms at 20 kHz
if ~isfield(o,'betaonesnap'), o.betaonesnap=10; end  % subsampling
if ~isfield(da,'timeseries'), da.timeseries=da.signal; da=rmfield(da,'signal'); warning('dataset signal field obsolete; use timeseries!'); end
if ~isfield(db,'timeseries'), db.timeseries=db.signal; db=rmfield(db,'signal'); warning('dataset signal field obsolete; use timeseries!'); end
if ~isfield(da,'name'), da.name='A'; end
if ~isfield(db,'name'), db.name='B'; end

Ya = arrayify(da.timeseries);  % sorting A...
Fa = arrayify(da.firings);
Ta = Fa(2,:); La = Fa(3,:); Ka = max(La);
Ca = ms_extract_clips2(Ya,Ta,o.T,[],o.betaonesnap);    % real t, resamples
Xa = ms_event_features(Ca,3);       % 3 features for viz
if o.verb, figure; set(gcf,'position',[1000 500 1000 1500]); tsubplot(3,2,1);
  ms_view_clusters(Xa,La); title([da.name ' labels in fea space']);
  Wopts.showcenter = 1;
  Wa = ms_templates(Ca,La);      % get mean waveforms (templates)
  tsubplot(3,2,2); ms_view_templates(Wa,Wopts); title([da.name ' templates']);
  drawnow
end

Yb = arrayify(db.timeseries);  % sorting B... (note timeseries could differ)
Fb = arrayify(db.firings);
if ~isempty(Fb)
  Tb = Fb(2,:); Lb = Fb(3,:); Kb = max(Lb);
  ii = Lb>0;        % B's classified labels
  Lb(~ii) = Kb+1; Kb = max(Lb);    % treat all unclass as one new label type
else, Tb = []; Lb = []; Kb = 0; end    % proceed with no-spike case
  
% match up B to A...
fprintf('best-permed accuracy of %s, treating %s as truth...\n',db.name,da.name)
if 1     % old matlab validspike way (3-pass)
  oo = o; oo.verb = 0;      % override output
  [perm Q acc t]=times_labels_accuracy(Ta,La,Tb,Lb,oo);
else     % C fast way, but doesn't do full times-labels best search
  outf = [tempdir '/accconfmat.mda'];    % lame for now
  mscmd_confusion_matrix(fnameify64(da.firings,tempdir),fnameify64(db.firings,tempdir),outf,o.max_matching_offset);   % todo: fix tempdir, don't overwrite
  Q = readmda(outf);   % un-permed
  [perm Q] = bestcolpermconfmat(Q);
  fprintf('best permed extended confusion matrix (%s down, %s across):\n',da.name, db.name)
  Q
  [~,acc] = labels_similarity(Q);
end
fk = acc.p;   % accuracy measure

popsa = histc(La,1:Ka);  % print info...
fprintf('%s populations for each label:\n',da.name);
fprintf('\t%d',1:Ka); fprintf('\n'); fprintf('\t%d',popsa); fprintf('\n');
popsb = histc(perm(Lb),1:Kb);
fprintf('%s populations for each label (best permuted):\n',db.name);
fprintf('\t%d',1:Kb); fprintf('\n'); fprintf('\t%d',popsb); fprintf('\n');
if isempty(Lb), warning('Lb labels are empty (no spikes found); no plots!');
  %if o.verb, close(gcf); end %removed by jfm 4/13/16
return; end

if o.verb
  Cb = ms_extract_clips2(Yb,Tb,o.T,[],o.betaonesnap);    % real t, resamples
  Xb = ms_event_features(Cb,3);       % 3 features for viz
  Wb = ms_templates(Cb,Lb);      % get mean waveforms (templates)
  % show B output now we know the best perm...
  tsubplot(3,2,3); ms_view_clusters(Xb,perm(Lb));
  title([db.name ' labels in fea space, best-permed']);
  [~,iperm] = sort(perm(1:Kb));          % invert permutation
  tsubplot(3,2,4); ms_view_templates(Wb(:,:,iperm),Wopts);
  title([db.name ' templates, best-permed']);

  % summarize confusion & accuracy...
  subplot(3,2,5); imagesc(Q); colorbar;ylabel([da.name ' label']);xlabel([db.name ' label']);
  hold on; plot([.5,Kb+.5;Kb+1.5,Kb+.5], [Ka+.5,.5;Ka+.5,Kb+1.5],'w-');
  title('best extended accuracy confusion matrix');
  subplot(3,2,6); plot(fk,'.','markersize',20); axis([1 Kb 0 1]);
  xlabel(sprintf('k (best permuted %s label)',db.name));
  ylabel('accuracy metric f_k');
end

if o.verb==2          % show timeseries and firings overlaid...
  addpath ~/spikespy/matlab/  % prefer old spikespy
  spikespy({Ya,Ta,La,sprintf('Y, %s',da.name)},{Yb,Tb,perm(Lb),sprintf('Y, %s',db.name)},{[t.tmiss';t.tfals';t.twrng'],[1+0*t.tmiss';2+0*t.tfals';3+0*t.twrng']});
  rmpath ~/spikespy/matlab/
end
if o.verb==3           % mountainview...
  mv.mode='overview2'; mv.raw=db.timeseries; mv.samplerate=db.samplerate;
  mv.firings=db.firings; mountainview(mv);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helpers...

function fname = fnameify32(X,outdir)
% FNAMEIFY  if array, writes to file and returns filename, otherwise keeps name

% v crude for now.

if nargin<2, outdir=tempdir; end
if ischar(X) || isstring(X)
  fname = X;
else
  fname = [outdir,'/',num2str(randi(1e10)),'.mda'];  % random filename
  writemda32(X,fname);
end

function fname = fnameify64(X,outdir)
% FNAMEIFY  if array, writes to file and returns filename, otherwise keeps name

% v crude for now.

if nargin<2, outdir=tempdir; end
if ischar(X) || isstring(X)
  fname = X;
else
  fname = [outdir,'/',num2str(randi(1e10)),'.mda'];  % random filename
  writemda64(X,fname);
end

function X = arrayify(X)
if ischar(X) || isstring(X), X = readmda(X); end


%%%%%%%%%%%
function test_compare_two_sortings   % tests default data and sorter
                                     % with two values of a parameter
d = demo_dataset;
o.samplerate = d.samplerate;
o.detect_threshold = 4.0;   % sorter opts: threshold in stddev EJ demo dataset
o.detect_polarity='m';
da.timeseries = d.timeseries;
da.outdir = [d.outdir 'a']; mkdir(da.outdir);
[da.firings,info] = simplesorter(da.timeseries,da.outdir,o);

o.detect_threshold = 3.5;    % change
db.timeseries = d.timeseries;
db.outdir = [d.outdir 'b']; mkdir(db.outdir);
[db.firings,info] = simplesorter(db.timeseries,db.outdir,o);

oo.verb = 0; fk = compare_two_sortings(da,db,oo);  % check different verbosities
oo.verb = 1; fk = compare_two_sortings(da,db,oo);

disp(['cleaning up temp dirs ' da.outdir ' and ' db.outdir ' ...'])
delete([da.outdir '/*']); rmdir(da.outdir);
delete([db.outdir '/*']); rmdir(db.outdir);
