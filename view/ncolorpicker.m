function c = ncolorpicker(n)
% NCOLORPICKER  n-by-3 RGB array of maximally perceptually distinct colors
%
% Barnett 6/9/16

if nargin<1, test_ncolorpicker; return; end

%c = ncolorpicker_repulsion(n);
c = ncolorpicker_snake(n);
%c = ncolorpicker_golden(n);

function c = ncolorpicker_golden(n)
% NCOLORPICKER  n-by-3 RGB maximally distinct colors by RGB hue golden jumps
% with property that sequence is the same for any n.
% Barnett 6/10/16

phi=(sqrt(5)-1)/2;
f = @(x) max(0,min(1,2-abs(x-3)));  % maps [0,6] -> [0,1]  for hue
t = 6*phi*(0:n-1)';
%a = mod(10000phi*(0:n-1),1)
%m = floor(mod(a*4,4));   % {0,1,2,3}
m = randi(4,n,1)-1;
% modulated speed: blue changes faster, and go thru RGB faster...
t = t + 0.2*sin(2*pi/3*(t-2.5)); % + 0.15*sin(pi*(t-1));
c = [f(mod(t,6)), f(mod(t-2,6)), f(mod(t+2,6))];   % snake on cube edges
%bri = rand(n,1).^.5; bri = 1+0*bri;
%sat = rand(n,1).^p; %sat = 1+0*sat;
%c = repmat(bri,[1 3]) .* (c.*repmat(sat,[1 3]) + repmat(1-sat,[1 3]));

for i=1:n                     % deterministic brightness and saturation cycle...
  if m(i)==1 || m(i)==3, c(i,:) = 0.65*c(i,:) + 0.35*[1 1 1]; end  % less sat for even
  if m(i)==1, c(i,:) = 0.7*c(i,:); end      % less brightness 4-cycle
end


function c = ncolorpicker_snake(n)
% NCOLORPICKER  n-by-3 RGB maximally distinct colors by RGB hue snake, cycle SV
%
% Barnett 6/9/16

%n = n-2; cgrey = [.5 .5 .5;1 1 1];   % kick off with grey & white, avoid from now
f = @(x) max(0,min(1,2-abs(x-3)));  % maps [0,6] -> [0,1]  for hue
t = 6*(1:n)'/n;
% modulated speed: blue changes faster, and go thru RGB faster...
t = 1.0 + t + 0.2*sin(2*pi/3*(t-2.5)) ; %+ 0.15*sin(pi*(t-1));
c = [f(mod(t,6)), f(mod(t-2,6)), f(mod(t+2,6))];   % snake on cube edges
%bri = rand(n,1).^.5; bri = 1+0*bri;
%sat = rand(n,1).^p; %sat = 1+0*sat;
%c = repmat(bri,[1 3]) .* (c.*repmat(sat,[1 3]) + repmat(1-sat,[1 3]));

for i=1:n                     % deterministic brightness and saturation cycle...
  if mod(i,2)==0, c(i,:) = 0.65*c(i,:) + 0.35*[1 1 1]; end  % less sat for even
  if mod(i,4)==0, c(i,:) = 0.7*c(i,:); end      % less brightness 4-cycle
end

%c = [c;cgrey]; n = n+2;  % append to grey stuff

% reorder so adjacent are far away...
%jump = max(primes(ceil(sqrt(n))));
if 1
phi=(sqrt(5)-1)/2; jump = max(primes(ceil(phi*n)));
i = 1+mod(jump*(0:n-1),n);
c = c(i,:);
end

% =============== OLD CODES AND THE SELF-TEST =================

function cbest = ncolorpicker_repulsion(n)
% NCOLORPICKER  n-by-3 RGB maximally distinct colors by Yukawa repulsion descent
%
% Barnett 6/9/16. Uses shorter metric in blue direction.

bluefac = 0.7;       % relative sensitivity to B vs R,G
cfixed = (0:0.1:0.3)'*ones(1,3);    % line of pts; black to grey to be avoided
cfixed(:,3) = cfixed(:,3)*bluefac;
restarts = 10;
niter = 50;
rate = 3e-3;
lam = 2.0 * n^(1/3);   % characteristic inverse length scale
Ebest = -inf;
for j=1:restarts
  c = rand(n,3);           % starting coords
  for i=1:niter
    [F E] = forcepotential([c;cfixed],lam);         % include fixed pts in calc
    %E, showc(c); %pause(0.01)
    c = c + rate*F(1:n,:);
    c(c>1) = 1; c(c<0) = 0;   % clip to unit cube
    c(c(:,3)>bluefac,3) = bluefac;   % clip blue        
    badinds = ~isfinite(sum(c,2)); c(badinds,:) = rand(sum(badinds),3);   % fix
    c(badinds,3) = c(badinds,3)*bluefac;
  end
  if E>Ebest, cbest = c; end
end
cbest(:,3) = cbest(:,3)/bluefac;
if ~isfinite(E), error('failed'); end

function [F V] = forcepotential(c,lambda)
% compute n-by-3 force vectors and scalar V, for Coulomb or Yukawa dist 1/lambda
n = size(c,1);
dx = repmat(c(:,1),[1 n]) - repmat(c(:,1)',[n 1]);
dy = repmat(c(:,2),[1 n]) - repmat(c(:,2)',[n 1]);
dz = repmat(c(:,3),[1 n]) - repmat(c(:,3)',[n 1]);
if 0               % Coulomb
  ir = 1./sqrt(dx.^2+dy.^2+dz.^2);
  ir(diagind(ir)) = 0;
  V = -sum(sum(triu(ir,1)));   % sum i>j -1/r interactions
  ir3 = ir.^3;
  F = [sum(dx.*ir3,2), sum(dy.*ir3,2), sum(dz.*ir3,2)];    % force rvec/r^3
else               % Yukawa
  r = sqrt(dx.^2+dy.^2+dz.^2); ir = 1./r;
  ir(diagind(ir)) = 0;
  e = exp(-lambda*r);
  V = -sum(sum(triu(ir.*e,1)));   % sum i>j -e/r interactions
  sc  = (1+lambda*r).*ir.^3.*e;    % scalar part
  F = [sum(dx.*sc,2), sum(dy.*sc,2), sum(dz.*sc,2)];    % force
end
  
function i = diagind(A)          % return indices of diagonal of square matrix
N = size(A,1); i = sub2ind(size(A), 1:N, 1:N);

function showc(c)
n = size(c,1);
hold off;
for j=1:n, plot3(c(j,1),c(j,2),c(j,3),'.','color',c(j,:),'markersize',20);
  hold on; text(c(j,1),c(j,2),c(j,3),sprintf('%d',j)); end
axis vis3d, axis([0 1 0 1 0 1]); drawnow

%%%%%%%%%%%%%%%
function test_ncolorpicker
n = 256;
c = ncolorpicker(n);
figure; showc(c)
figure; image(reshape(c,[1,n,3]));
% output for C static:

fprintf('float c[256][3] = {\n')
for i=1:n
  fprintf('{ %.3f, %.3f, %.3f },\n',c(i,1),c(i,2),c(i,3))
end
fprintf('};\n')
