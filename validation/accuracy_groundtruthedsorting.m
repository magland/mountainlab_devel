function [fk Q perm] = accuracy_groundtruthedsorting(dataset,o_acc)
% ACCURACY_GROUNDTRUTHEDSORTING  meas accuracy of sorter on gnd-truth data
%
% fk = accuracy_anysorter_groundtrutheddata(sortfunc,dataset,o_acc,o_sorter)
%   measures any sorter's accuracy vs known firing times on synthetic demo data.
%   Also pops up a bunch of windows showing accuracy.
% [fk Q perm] = ... also outputs comparison info.
%
% Inputs:
%  sortfunc - function handle to a sorter with the interface
%                [firingsfile,info] = sorter(rawfile,outdir);
%  dataset - struct with fields giving dataset with its ground truth:
%            dataset.timeseries - MDA filename or M*Ns array, raw/filtered timeseries
%            dataset.firings - MDA filename or 4*Ns array, sorted firings
%            dataset.truefirings - MDA filename or 4*Ns array, true firings
%            dataset.truewaveforms - MDA filename or M*T*K array, true waveforms
%                                    (optional)
%            dataset.samplerate - samplerate in samples/sec for raw data
%            dataset.name - (optional) string name of dataset
%  o_acc - accuracy measuring options:
%            o_acc.usepre          - if true, use preprocessed (pre.mda) not raw
%                                    (default false)
%            o_acc.xc - if true, show cross-correlograms
%            (others passed to compare_two_sortings; see its help)
% Outputs:
%  fk - (1xK) accuracy metrics on the neuron types labeled by ground truth labels
%  Q, perm - extended confusion matrix and best perm, as output by
%            compare_two_sortings
%
% Notes: mostly a wrapper to compare_two_sortings

% Barnett 7/21/16 based on accuracy_anysorter_groundtrutheddata

if nargin<2, o_acc=[]; end
if ~isfield(dataset,'name'), dataset.name = ''; end

% note timeseries is ignored
da = dataset; db = da;   % now just what needs to change...
da.name = [dataset.name ' true'];
da.firings = dataset.truefirings;
db.name = [dataset.name ' sorted'];            % dataset B is the sorting

[fk Q perm] = compare_two_sortings(da,db,o_acc);

if isfield(o_acc,'xc') && o_acc.xc % do cross-corr plot....
  f = readmda(db.firings);   % arrayify?
  t = f(2,:); l = perm(f(3,:));       % so labels match those in compare plots
  show_crosscorr(l,t);  % is v slow to plot
end
