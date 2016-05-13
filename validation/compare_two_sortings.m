function [fk Q perm iperm] = compare_two_sortings(da,db,o)
% COMPARE_TWO_SORTINGS  accuracy metrics between two sortings treating one as true
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
%      ts - (default true). If true, read timeseries, otherwise just firings
%      rowcol - (default true), if verb>0 output 2nd plot of row- & col-normed Q
%      T - clip size in samples
%      max_matching_offset - passed to accuracy measuring
%      betaonesnap - passed to clip extraction
%      permfirings (default true) - output permuted B firings
%      verb - verbosity: 0 (text only), 1 (figure), 2 (fig+spikespy),
%                        3 (fig+MV), 4 (fig+clipclouds).
% Output:
%  fk - (1xK) accuracy metrics of d2 treating d1 as ground truth.
%  Q - extended best-permuted confusion matrix
%  perm - best permutation of B's labels to match A, ie LA ~ perm(LB)
%         (if there are no types k in LB, kth entry of perm is NaN).
%  iperm - best permutation of A's labels to match B.
%
% Without arguments, does self-test
%
% Used by: accuracy_anysorter_groundtrutheddata.m

% todo: opt to switch to mscmd_conf_mat for speed, test.
% Barnett 4/5/16. Kb absent label issue 4/19/16. opts.ts 5/12/16

if nargin==0, test_compare_two_sortings; return; end
if nargin<3, o = []; end
if ~isfield(o,'verb'), o.verb = 1; end
if ~isfield(o,'ts'), o.ts = 1; end
if ~isfield(o,'rowcol'), o.rowcol = 1; end
if ~isfield(o,'permfirings'), o.permfirings = 1; end
if ~isfield(o,'T'), o.T = 50; end
if ~isfield(o,'max_matching_offset'), o.max_matching_offset = 10; end  % in samples; eg 0.5 ms at 20 kHz
if ~isfield(o,'betaonesnap'), o.betaonesnap=10; end  % subsampling
if ~isfield(da,'timeseries'), da.timeseries=da.signal; da=rmfield(da,'signal'); warning('dataset signal field obsolete; use timeseries!'); end
if ~isfield(db,'timeseries'), db.timeseries=db.signal; db=rmfield(db,'signal'); warning('dataset signal field obsolete; use timeseries!'); end
if ~isfield(da,'name'), da.name='A'; end
if ~isfield(db,'name'), db.name='B'; end

Fa = arrayify(da.firings);                      % sorting A...
Ta = Fa(2,:); La = Fa(3,:); Ka = max(La);
if o.ts, Ya = arrayify(da.timeseries); 
  Ca = ms_extract_clips2(Ya,Ta,o.T,[],o.betaonesnap);    % real t, resamples
  Xa = ms_event_features(Ca,3);       % 3 features for viz
end
if o.verb, bigfig=figure;
  if o.ts, set(gcf,'position',[1000 500 1000 1500]);
    tsubplot(3,2,1);
    ms_view_clusters(Xa,La); title([da.name ' labels in fea space']);
    Wopts.showcenter = 1; Wopts.pops = histc(La,1:max(La));  % cluster sizes
    Wa = ms_templates(Ca,La);      % get mean waveforms (templates)
    tsubplot(3,2,2); ms_view_templates(Wa,Wopts); title([da.name ' templates']);
    drawnow
  else, set(gcf,'position',[1000 500 2000 1000]);
  end
end

Fb = arrayify(db.firings);    % sorting B... (note timeseries could differ)
if ~isempty(Fb)
  Tb = Fb(2,:); Lb = Fb(3,:); Kb = max(Lb); % K = # label types in sorting
  ii = Lb>0;        % B's classified labels
  Lb(~ii) = Kb+1; Kb = max(Lb);    % treat all unclass as one new label type
else, Tb = []; Lb = []; Kb = 0; end    % proceed with no-spike case

% match up B to A...
fprintf('best-permed accuracy of %s, treating %s as truth...\n',db.name,da.name)
if 1     % old matlab validspike way (3-pass)
  oo = o; oo.verb = 0;      % override output
  [perm Q acc t]=times_labels_accuracy(Ta,La,Tb,Lb,oo);
else     % C fast way, but doesn't do full times-labels best search. todo: use
  outf = [tempdir '/accconfmat.mda'];    % lame for now
  mscmd_confusion_matrix(fnameify64(da.firings,tempdir),fnameify64(db.firings,tempdir),outf,o.max_matching_offset);   % todo: fix tempdir, don't overwrite
  Q = readmda(outf);   % un-permed
  [perm Q] = bestcolpermconfmat(Q);
  fprintf('best permed extended confusion matrix (%s down, %s across):\n',da.name, db.name)
  Q
  [~,acc] = labels_similarity(Q);
end
fk = acc.p;   % accuracy measure
[kblist,iperm] = sort(perm(1:Kb));  % best-permed labels present; inverse perm
meaningful = ~isnan(kblist);
kblist = kblist(meaningful); iperm = iperm(meaningful);
mKb = max(kblist);   % largest permed B type, should match size(Q,2)-1

popsa = histc(La,1:Ka);  % print info...
fprintf('%s populations for each label:\n',da.name);
fprintf('\t%d',1:Ka); fprintf('\n'); fprintf('\t%d',popsa); fprintf('\n');
popsbwempties = histc(perm(Lb),1:mKb);   % B's permed pops including empty labels
popsb = popsbwempties(kblist);   % B's permed pops only in the kblist labels
fprintf('%s populations for each label (best permuted):\n',db.name);
fprintf('\t%d',kblist); fprintf('\n'); fprintf('\t%d',popsb); fprintf('\n');
if isempty(Lb), warning('Lb labels are empty (no spikes found)');
  if o.verb, close(bigfig); end % removed by jfm 4/13/16; close only g, ahb 4/19
return; end
if o.permfirings
  [path,nam] = fileparts(db.firings);
  pffile = [path,'/',nam,'_permed.mda'];
  try
    writemda64([0*Tb;Tb;perm(Lb)],pffile);  % output best-permed B firings
  catch
    warning('writemda failed');
  end
end

if o.verb
  if o.ts, Yb = arrayify(db.timeseries);
    Cb = ms_extract_clips2(Yb,Tb,o.T,[],o.betaonesnap);    % real t, resamples
    Xb = ms_event_features(Cb,3);       % 3 features for viz
    Wb = ms_templates(Cb,perm(Lb));     % get mean waveforms (templates) up to mKb
    % show B output now we know the best perm...
    tsubplot(3,2,3); ms_view_clusters(Xb,perm(Lb));
    title([db.name ' labels in fea space, best-permed']);
    Wopts.pops = popsbwempties;
    tsubplot(3,2,4); ms_view_templates(Wb,Wopts);  % NB Wb already permed
    title([db.name ' templates, best-permed']);
  end

  % summarize confusion & accuracy...
  if o.ts, subplot(3,2,5); else, subplot(1,2,1); end
  imagesc(Q); colorbar; title('best extended accuracy confusion matrix');
  ylabel([da.name ' label']);xlabel([db.name ' label']);
  if o.ts, subplot(3,2,6); else, subplot(1,2,2); end
  bar(fk); axis([0.5 mKb+0.5 0 1]);
  xlabel(sprintf('k (best permuted %s label)',db.name));
  ylabel('accuracy metric f_k');
end
drawnow

if o.verb & o.rowcol      % 2nd figure w/ more detail on normed Q etc.....
  display_confusion(Q,da.name,db.name);
end

if o.ts %
if o.verb==2          % show timeseries and firings overlaid...
  %addpath ~/spikespy/matlab/  % prefer old spikespy
  spikespy({Ya,Ta,La,sprintf('Y, %s',da.name)},{Yb,Tb,perm(Lb),sprintf('Y, %s',db.name)},{[t.tmiss';t.tfals';t.twrng'],[1+0*t.tmiss';2+0*t.tfals';3+0*t.twrng']});
  %rmpath ~/spikespy/matlab/
end
if o.verb==3           % mountainview... (not best-permuted!)
  mv.mode='overview2'; mv.raw=db.timeseries; mv.samplerate=db.samplerate;
  mv.firings=db.firings; mountainview(mv);
end
if o.verb==4           % clip clouds...
  cco=[]; cco.vzoom = 1.0;
  figure; ms_view_clip_clouds(Ca,La,cco); set(gcf,'name','true clips');
  figure; ms_view_clip_clouds(Cb,perm(Lb),cco); set(gcf,'name','sorted clips');
  figure(bigfig);
end
end %

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

%oo.verb = 0; fk = compare_two_sortings(da,db,oo);  % check different verbosities
%oo.verb = 1; fk = compare_two_sortings(da,db,oo);
oo.verb = 1; oo.ts = 0; [fk Q perm iperm] = compare_two_sortings(da,db,oo);
perm
ls(db.outdir)  % check firings_permed.mat is there

disp(['cleaning up temp dirs ' da.outdir ' and ' db.outdir ' ...'])
delete([da.outdir '/*']); rmdir(da.outdir);
delete([db.outdir '/*']); rmdir(db.outdir);
