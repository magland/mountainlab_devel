function ss_view_clusters(z,L,opts)
% Snapshot of plot_labeled_points on 5/22/2015
% SS_VIEW_CLUSTERS - shows points in 2d or 3d space colored by their labels
%
% ss_view_clusters(z,L) where z is 2-by-N or 3-by-N and L is 1-by-N
%  with entries in 1...K, plots each col of z as a point in R^2 or R^3 as
%  appropriate, with color given by corresponding entry in L. If z has >3 rows,
%  just the first 3 are used.
%  Unclassified pits (L=0) are shows as grey plus signs.
%
% Barnett 12/18/14. Added grey "+" 12/19/14

addpath([fileparts(mfilename('fullpath')),'/internal']);
addpath([fileparts(mfilename('fullpath')),'/internal/colorspace']);

warning('off','MATLAB:legend:IgnoringExtraEntries'); %by jfm 2/12/15 -- warning from legnum at bottom

if nargin<3 opts=struct; end;
if (~isfield(opts,'marker_size')) opts.marker_size=5; end;
if (~isfield(opts,'create_figure')) opts.create_figure=1; end;

K = max(L);
[dims N] = size(z);
if N~=numel(L), error('length of L must match number of columns of z!'); end
%c = 'bgrkmcy';  % colors repeat if >7 labels
if (opts.create_figure) figure; end;

CC=distinguishable_colors(K,{'w'});
colors={};
for j=1:size(CC,1)
	colors{j}=CC(j,:);
end;


marker_size=opts.marker_size;

if dims>=3
  for n=1:K, j=L==n; % for some obscure reason an empty plot3 locks to 2d view!
    if sum(j), plot3(z(1,j),z(2,j),z(3,j),'.','color',colors{n}, 'markersize',marker_size);
      hold on; end % show labels via color
  end
  %j=L==0; plot3(z(1,j),z(2,j),z(3,j), '+', 'color',[.5 .5 .5],'MarkerSize',marker_size);
  axis equal vis3d; xlabel('z_1');ylabel('z_2');zlabel('z_3');

elseif dims==2
  for n=1:K, j=L==n; plot(z(1,j),z(2,j),[c(mod(n-1,numel(c))+1) '.']);
    hold on; end % show labels via color
    j=L==0; plot(z(1,j),z(2,j), '+', 'color',[.5 .5 .5]);
  axis equal; xlabel('z_1');ylabel('z_2');
  
else, error('size(z,1) must be 2 or more!');
end

hline(0); vline(0); % makes origin obvious
legnum(1:K);     % label colors

hold off; %added by jfm 10/29/15


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

