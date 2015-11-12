function ss_view_waveforms(w,title0,vertical_spread,t,fig,o)
% Snapshot of plot_spike_shapes on 5/22/2015
% ss_view_waveforms - plot a set of multi-channel waveforms or clips
%
% ss_view_waveforms(W) where W is a M x Nt x Ns array with M channels, Nt
%  time samples, and Ns spikes (or events) plots a 2D grid of graphs,
%  with the channel increasing vertically downwards, and the spike number
%  increasing horizontally.
%
% jfm added optional fig parameter (handle to figure to use) - 2/11/15
%
% plot_spike_shapes(W, title)
% plot_spike_shapes(W, title, vertical_spread)
%  Controls the title for the plot, and the vertical separation between
%  graphs in the signal units.
%
% plot_spike_shapes(W, title, vertical_spread, ta) adds alignment lines at
%  the times given by ta, measured in offset from the first time index.
%  If ta is a numeric array, one line per clip is shown.
%  If ta is a struct array, the times ta(j).t are shown in clip j (ie ta is
%   treated as a parameter array).
%
% If W is a struct, it's treated as variable-length clip object with fields
%  X, Ts, tptr, etc (see mergeclips).
%
% plot_spike_shapes(W, title, vertical_spread, ta, fig) uses existing figure
%  number.
% 
% plot_spike_shapes(W, title, vertical_spread, ta, [], opts) control options:
%    opts.lines = 0,1 switches veritcal dotted lines from t.t
%
% todo: If the field d is present in W, a 1 ms time bar added
%
% Without input arguments it runs a self-test using data in data/

% todo: various tweaks such as allowing axes, figure placement?
% todo: index set input to show subset w/ correct clip # labels.
%
% Magland/Barnett
% 12/16/14 ahb: added self-test, tidied up. 1/23/15 annotation, autosize
% 1/29/15 added optional alignment time plotting. 2/19/15: t may be struct array
% 2/19/15 simplified padding mess, and W may be struct for variable-width clips
% 2/27/15 opts, lines, existfig
% 3/31/15 vertical_spread can be <0 for flipping ordering

if nargin<6, o = []; end               % opts
if ~isfield(o,'lines'), o.lines = 1; end
if ~isfield(o,'create_figure'), o.create_figure = 1; end

if (nargin>=2)
	if (isstruct(title0)) o=title0; title0=[]; end;
end;

w=w(end:-1:1,:,:); % Let's re-order things so we can stay sane! jfm 6/4/2015

if nargin<1, test_view_waveforms; return; end
if nargin<2 || isempty(title0), title0=''; end
if ~isstruct(w), w = mergeclips(w); end   % convert to struct in variable-len
if nargin<3 || isempty(vertical_spread), vertical_spread = 1.0 * max(abs(w.X(:))); end; % auto scale
if vertical_spread==0, vertical_spread = 1e-16; end % tiny value (in case w.X=0)
if nargin<4, t=[]; end
if (~o.create_figure)
	existfig=1; fh=gcf;
else
	if (nargin<5 || isempty(fig)) existfig=0; fh=figure;
	else existfig=1; fh=fig; figure(fh); end
end;

padding=2;           % means that 1st sample starts at x-coord = padding+1 (=3)
w = padclips(w,padding,nan);
NC = size(w.X,1); % spread out the channels
for c=1:NC
	w.X(c,:)=w.X(c,:)+c*vertical_spread;
end
sx=min(3000,50+0.7*w.Ttot); sy=320;                % controls figure size

plot(w.X'); title(title0);        % show data

if ~existfig, set(fh,'Position', [100,100,sx,sy]); end
xlim([1,w.Ttot]);
ylim(sort(vertical_spread*[.3,NC+.8+0.05*NC])); % handles either sign

K = w.Ns; % Alex's annotations:
if K>1   % number the clips
  for k=1:K, text(w.tptr(k) + w.Ts(k)/2, vertical_spread*(NC+0.6), sprintf('%d',k), 'rotation',90,'color',.7*[1 1 1]); end
end
if (length(w.tptr)>1) %condition added by jfm (2/19/15)
	h=vline(w.tptr(2:end), '-'); set(h,'color',.85*[1 1 1]);   % faint dividers
end;
if ~isempty(t)                  % show alignment times
  for k=1:K, p = t(k);          % get param struct or scalar
    if ~isstruct(p), p.t = p; end  % fake a param struct
    annotateparams(p, w.tptr(k)+padding, -0.2*vertical_spread, o);
  end
end
set(fh,'Color',[1,1,1]);
set(gca,'YTickLabel',[]); set(gca,'Ytick',[]);  % remove axes
set(gca,'XTickLabel',[]); set(gca,'Xtick',[]);
% jfm: where is box set on?
%%%%%%

function annotateparams(p,ioff,ytxt,o)
% show parameters p.t, p.l, offset by index ioff in x. ytxt gives text y
if o.lines, h=vline(ioff + p.t, '--'); set(h,'color',.5*[1 1 1]); end
if isfield(p,'l')
  for j=1:numel(p.l), h=text(ioff + p.t(j), ytxt, sprintf('%d',p.l(j))); end
end


function m = mergeclips(A,B)
% MERGECLIPS - merge two sets of variable- or fixed-length clip arrays X
%
% m = mergeclips(A,B) produces a variable-length clip struct m from A & B,
%  in that order.
%  A, B may each be either 3d (fixed-length) clips or variable-length clip
%  structs.
%
% The fixed-length format is
%    X (double M*Nt*Ns) - signal array. M channels, Nt time points
%        per clip, and Ns clips.
%
% The variable-length clip struct m (or inputs A and/or B) has fields:
%    X (double M*Ttot) - concatenated signal array. M channels
%    Ttot (int) - total time points = sum(Ts)
%    Ts (int 1*Ns) - clip lengths in time points
%    tptr (int 1*Ns) - indices of starts of each clip
%    Ns (int) - number of clips
%
% m = mergeclips(X) produces variable-length clip struct from 3D array X.

% todo: think about if A, B huge in RAM.
%
% Barnett 2/19/15

if nargin<2, B.X = []; B.Ts = []; end  % dummy B so A converted to m

if isstruct(A) && isstruct(B)
  m.X = [A.X B.X]; m.Ts = [A.Ts B.Ts];
elseif isstruct(A) && ~isstruct(B)
  [M Nt Ns] = size(B);
  m.X = [A.X reshape(B,[M Nt*Ns])];
  m.Ts = [A.Ts Nt*ones(1,Ns)];
elseif ~isstruct(A) && isstruct(B)
  [M Nt Ns] = size(A);
  m.X = [reshape(A,[M Nt*Ns]), B.X];
  m.Ts = [Nt*ones(1,Ns) B.Ts];
else
  [M NtA NsA] = size(A); [M NtB NsB] = size(B);
  m.X = [A B]; m.Ts = [NtA*ones(1,NsA) NtB*ones(1,NsB)];
end
m.Ttot = sum(m.Ts);                 % less to go wrong than if case-by-case
m.Ns = numel(m.Ts);
m.tptr = [1 1+cumsum(m.Ts(1:m.Ns-1))];

function n = padclips(m,padding,padval)
% PADCLIPS - pad along the time axis a struct of variable-length clips
%
% B = padclips(A,padding,padval) pads with padval the data A.X in A to give B.X,
%  using "padding" on each end of the time axis of each clip.
%
% See also: plot_spike_shapes, which uses it; mergeclips for struct format

% Barnett 2/19/15
n = m;   % copies over unknown fields, and Ns
n.Ts = m.Ts + 2*padding;
n.Ttot = sum(n.Ts);
n.tptr = [1 1+cumsum(n.Ts(1:n.Ns-1))];
M = size(m.X,1);
n.X = padval*ones(M,n.Ttot);   % start a fresh n.X. JFM's tmp array
for c=1:m.Ns               % could do in one go like subsetclips
  j = 0:m.Ts(c)-1;
  n.X(:,n.tptr(c)+padding+j) = m.X(:,m.tptr(c)+j);  
end


%%%%%%%
function test_view_waveforms   % needs waveforms data in data/
for c='hbe'
  load(sprintf('data/waveforms_%s_fac1',c));  % get not upsampled data
  [M NT K] = size(W); t = ones(1,K) * (NT-1)/2;  % std t_aligns for wf
  plot_spike_shapes(W,c,[],t);  % alignments should be visually correct
end
% todo: test variable-length clips