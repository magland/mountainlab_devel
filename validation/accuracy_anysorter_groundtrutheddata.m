function [fk Q perm info] = accuracy_anysorter_groundtrutheddata(sortfunc,dataset,o_acc,o_sorter)
% ACCURACY_ANYSORTER_GROUNDTRUTHEDDATA  meas accuracy of sorter on gnd-truth data
%
% fk = accuracy_anysorter_groundtrutheddata(sortfunc,dataset,o_acc,o_sorter)
%   measures any sorter's accuracy vs known firing times on synthetic demo data.
%   Also pops up a bunch of windows showing accuracy.
% [fk Q perm info] = ... also outputs comparison info.
%
% Inputs:
%  sortfunc - function handle to a sorter with the interface
%                [firingsfile,info] = sorter(rawfile,outdir);
%  dataset - struct with fields giving dataset with its ground truth:
%            dataset.timeseries - MDA filename or M*N array, raw EC signal data
%            dataset.truefirings - MDA filename or 4*Ns array, true firings
%            dataset.truewaveforms - MDA filename or M*T*K array, true waveforms
%                                    (optional)
%            dataset.samplerate - samplerate in samples/sec for raw data
%            dataset.outdir - path to directory for MDA output (needn't exist)
%            dataset.name - (optional) string name of dataset
%  o_acc - accuracy measuring options:
%            o_acc.usepre          - if true, use preprocessed (pre.mda) not raw
%                                    (default false)
%            o_acc.xc - if true, show cross-correlograms
%            (others passed to compare_two_sortings; see its help)
%  o_sorter (optional) - passed to last argument of sortfunc, as in:
%                [firingsfile,info] = sortfunc(rawfile,outdir,o_sorter)
%                Note: the sorter opts will already be passed the samplerate
%                      as given in the dataset struct.
% Outputs:
%  fk - (1xK) accuracy metrics on the neuron types labeled by ground truth labels
%  Q, perm - extended confusion matrix and best perm, as output by
%            compare_two_sortings
%  info -    any info struct returned by the sorter.
%
% Run without options runs a test on demo data (accuracy_simplesorter)

% Barnett 3/3/16 based on accuracy_sort_demotimeseries.m
% 3/16/16 Barnett changed to use non-integer firing time resampling
% todo: * make various switches below into options.
% 3/18/16 struct not text input for dataset. truewaveforms optional
% 4/6/16 vastly simplified by using compare_two_sortings. 4/8/16 prefile

if nargin==0, accuracy_simplesorter; return; end             % is in validation/
if nargin<2|isempty(dataset), dataset = demo_dataset; dataset.name = 'demo'; end
if nargin<3, o_acc=[]; end
if nargin<4, o_sorter=[]; end
if ~isfield(dataset,'name'), dataset.name = ''; end

outdir = dataset.outdir; if ~exist(outdir,'dir'), mkdir(outdir); end
samplerate = dataset.samplerate;
o_sorter.samplerate = samplerate;

da.timeseries = fnameify32(dataset.timeseries,outdir);  % dataset A struct "true"
da.firings = dataset.truefirings; da.name = [dataset.name ' true'];
db = da; db.name = [dataset.name ' sorted'];            % dataset B, from sorter
db.samplerate = samplerate;
% note da and db just point to the same timeseries but different firings files

disp('call sorter:')
[db.firings,info] = sortfunc(da.timeseries,outdir,o_sorter);

%mv.pre = info.prefile; mv.firings=db.firings; mountainview(mv); return
%%%%%%%%%%%%%%%%%%% debugging

if isfield(o_acc,'usepre') & o_acc.usepre
  if exist(info.prefile,'file')
    disp('using pre-proc (not raw) timeseries for comparsion plots...');
    da.timeseries = info.prefile; db.timeseries = da.timeseries;
  else
    warning('no pre-proc file found, using raw timeseries!');
  end
end
  
[fk Q perm] = compare_two_sortings(da,db,o_acc);
%disp('accuracies:'), fk

if isfield(o_acc,'xc') && o_acc.xc % do cross-corr plot....
  f = readmda(db.firings);   % arrayify?
  t = f(2,:); l = f(3,:);
  show_crosscorr(l,t);  % is v slow to plot
end

%%%%%%%%%%%%%%%

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
%%%%%%%%
