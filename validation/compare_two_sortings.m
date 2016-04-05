function fk = compare_two_sortings(d1,d2,o)
% COMPARE_TWO_SORTINGS  accuracy metrics between two sortings treating one true
%
% fk = compare_two_sortings(d1,d2,o)
% produces several figures and text output.
%
% Inputs:
%  d1, d2 - each a dataset struct with at least the following fields:
%           timeseries - MDA filename or M*N array, raw EC signal
%           firings    - MDA filename or 4*Ns array, firings
%                        (row 1 is peak channels, row 2 firing times t_j, row 3
%                        is firing identities k_j, row 4 is firing amplitudes)
% Output:
%  fk - (1xK) accuracy metrics of d2 treating d1 as ground truth.
%
% To be used by: accuracy_anysorter_groundtrutheddata.m
%
% Barnett 4/5/16

if ~isfield(d1,'timeseries'), d1.timeseries=d1.signal; d1=rmfield(d1,'signal'); warning('dataset signal field obsolete; use timeseries!'); end
if ~isfield(d2,'timeseries'), d2.timeseries=d2.signal; d2=rmfield(d2,'signal'); warning('dataset signal field obsolete; use timeseries!'); end

Y1 = arrayify(d1.timeseries);
F1 = arrayify(d1.firings);
Y2 = arrayify(d2.timeseries);
F2 = arrayify(d2.firings);

% ...











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
