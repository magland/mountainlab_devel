function ms_view_clusters(features,labels)
%MS_VIEW_CLUSTERS - 2D or 3D view of clusters in feature space, colored by
%labels
%
% Syntax:  ms_view_clusters(features,labels)
%
% Inputs:
%    features - RxL array of features (R>=2), e.g., optained by
%               mscmd_features or ms_event_features
%    labels - 1xL array of integer labels (controls the colors of the data
%    points)
%
% Other m-files required: distinguishable_colors
%
% See also: mscmd_features, ms_event_features

% Author: Jeremy Magland
% Jan 2015; Last revision: 15-Feb-2106

if nargin<1, test_ms_view_clusters; return; end;

addpath([fileparts(mfilename('fullpath')),'/colorspace']);

M=size(features,1); if (M>3) features=features(1:3,:); end; M=size(features,1);
N=size(features,2);

if nargin<2
    labels=ones(1,N);
end;

K=max(labels);

CC=distinguishable_colors(K,{'w'});
colors={};
for j=1:size(CC,1)
	colors{j}=CC(j,:);
end;

if M==2
    for k=0:K
        inds=find(labels==k);
        if (length(inds)>0)
            if (k>0)
                plot(features(1,inds),features(2,inds),'.','Color',colors{k}); hold on;
            else
                plot(features(1,inds),features(2,inds),'.','Color',[0.5,0.5,0.5]); hold on;
            end;
        end;
    end;
    hold off;
    legnum(1:K);
elseif M==3
    for k=0:K
        inds=find(labels==k);
        if (length(inds)>0)
            if (k>0)
                plot3(features(1,inds),features(2,inds),features(3,inds),'.','Color',colors{k}); hold on;
            else
                plot3(features(1,inds),features(2,inds),features(3,inds),'+','Color',[0.5,0.5,0.5]); hold on;
            end;
        end;
    end;
    hold off;
    legnum(1:K);
else
    error('Invalid number of dimensions: %d',M);
end

end

function test_ms_view_clusters
X=cat(2,randn(2,200),3+randn(2,200));
labels=cat(2,ones(1,200),ones(1,200)*2);
figure;
ms_view_clusters(X,labels);
end

function legnum(a, prec, prefix)

if nargin==1
  legend(num2cellstr(a));
elseif nargin==2
  legend(num2cellstr(a, prec));
elseif nargin==3
  legend(num2cellstr(a, prec, prefix));
else
  error('too many arguments to legnum.')
end
end

%
% NUM2CELLSTR convert array of floating-point numbers to cell array of strings
%    NUM2CELLSTR(X) converts array X to cell array of strings.
%    If X is a two- or multi-dimensional array, it will be
%    flattened (all elements will still be included).
%
%    NUM2CELLSTR(X, P) is the same but uses precision P, where P is an integer.
%
%    NUM2CELLSTR(X, P, S) same as above but includes a prefix string to
%    each cell.
%
%    This clumsy routine would be unnecessary if Matlab provided something
%    like python's string.strip() function.
%
% See also SPRINTF, CELLSTR
%
%    Alex Barnett 12/5/02

function [c] = num2cellstr(a, prec, prefix)

if nargin==1
  prec = 4;      % default precision
else
  if prec<1
    error('precision must be at least 1.')
  end
  if prec>16
    error('precision cannot exceed 16.')
  end
end
if nargin<3
  prefix = ''; % default prefix
end

l = 25;          % max number of characters for representing a number
n = numel(a);

% build printf format string
f = sprintf('%%-%d.%dg', l, round(prec));

c = cellstr([repmat(prefix, [n 1]) reshape(sprintf(f, a),[l, n])']);
end


%
% LEGNUM Legend current figure using array of numbers.
%    LEGNUM(X) adds a legend to current figure using string
%    representations of the numbers in X. If X is a two- or multi-dimensional
%    array, it will be flattened and all elements will be included.
%
%    LEGNUM(X, P) is the same but uses precision P, where P is an integer.
%
%    LEGNUM(X, P, S) same as above but includes a prefix string to
%    each legend label.
%
% Examples
%    legnum(logspace(-5,-4,7), 6);
%    Adds a legend with logarithmically-spaced number labels, with
%    6 significant digit precision
%
%    legnum(logspace(-5,-4,7), 6, 'x = ');
%    Same but labels are of the form 'x = 1e-5', etc.
%
% See also NUM2CELLSTR
%
%    Alex Barnett 12/5/02

function legnum(a, prec, prefix)

if nargin==1
  legend(num2cellstr(a));
elseif nargin==2
  legend(num2cellstr(a, prec));
elseif nargin==3
  legend(num2cellstr(a, prec, prefix));
else
  error('too many arguments to legnum.')
end
end