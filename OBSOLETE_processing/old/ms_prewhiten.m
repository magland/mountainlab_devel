function X=ms_prewhiten(X,opts)
X=channelprewhiten(X,opts.samplefreq);
end

function A=channelprewhiten(A,samplefreq)

% channelprewhiten - spatially pre-whiten a timeseries using C from noise clips
%
% function A = channelprewhiten(A,thresh,o)
%  A is timeseries data, samplefreq is the sampling frequency in Hz
%  On output A is changed.

C = empiricalspacetimecorr(A,samplefreq,[]);
Q=inv(chol(squeeze(C(:,:,(size(C,3)+1)/2)))');   % matrix to transform so
% rows of d.A are indep and unit-norm
for i=1:size(Q,1), Q(i,:) = Q(i,:)/norm(Q(i,:)); end  % preserve norms of rows
A=Q*A;     % transform data

end

function C = empiricalspacetimecorr(A,samplefreq,thresh,o);
% EMPIRICALSPACETIMECORR  estimate "noise" autocorrelation over channels and time
%
% C = empiricalspacetimecorr(d,thresh,o)
%
% Inputs:
%  A -M-by-M data, samplefreq sampling time
%  thresh - (optional) threshold below which counted as noise clip, measured if
%           not given or empty.
% Outputs:
%  C - M*M*Nt autocorrelation matrix
%
% dtau is fixed at 1 time sample for now

dt=1/samplefreq;

if nargin<3 || isempty(thresh), thresh = autothreshold(A,samplefreq); end
if nargin<4, o=struct; end
if ~isfield(o,'Twin'), o.Twin = 0.004; end % window length for below-thresh clips
%if ~isfield(o,'Twin'), o.Twin = 0; end % window length for below-thresh clips
if ~isfield(o,'verb'), o.verb = 0; end
Nt = ceil(o.Twin/dt); if mod(Nt,2)==0, Nt=Nt+1; end % odd, window width
maxtau = (Nt-1)/2; taus = -maxtau:maxtau;  % dtau fixed at 1 sample
[M N] = size(A);

fprintf('computing space-time cross-corr... ')
n = 1e4;  % # noise-clip trials
i=1; C = zeros(M,M,Nt);  % cross-corr
while i<=n
  j = randi(N,1);
  if j>maxtau && j<=N-maxtau
    X = A(:,j+taus);
    if min(X(:))>-thresh   % it's noise, not an "event"
      for m=1:M
        C(m,:,:) = permute(C(m,:,:),[2,3,1]) + X(m,1+maxtau)*X; % 1+maxtau is central el, jfm replaced squeeze by permute to handle case of Nt=1
      end
      i = i+1;
      if mod(i,round(n/10))==0, fprintf('%d%% ',round(100*i/n)); end
    end
  end  
end
C = C/n;
fprintf('\n')

end

function t = autothreshold(A,samplefreq)
% AUTOTHRESHOLD  choose threshold for detection given raw EC data struct d
%  containing d.A timeseries, d.samplefreq
%
% function t = autothreshold(d) returns a single threshold based on a heuristic
%  noise model

% Barnett 6/11/15

noi = empiricalnoise(A,samplefreq);
%fprintf('\tnoi.eta = %.3g\n',noi.eta)
t = 5.0*noi.eta;   % how many std deviations from zero
end


function noi = empiricalnoise(A,samplefreq)
% EMPIRICALNOISE - extract noise model parameters from raw EC dataset
%
% noi = empiricalnoise(A). Currently two methods for noi.eta built in.
%
% Inputs:
%  d - raw EC dataset struct, with fields: A - (M*Nt) raw data
%                                          dt - timestep
%                                          etc
% Outputs:
%  noi - noise model struct with fields: eta - std error in iid Gauss model
%
% todo: self-test

% Barnett 2/12/15 

meth='j';   % choose method.   todo: make opt

if strcmp(meth,'a')          % Alex idea: curvature of small signal data
  sc = max(abs(A(:))); % fit a Gaussian noise model to small-ampl data...
  b = (-1:0.01:1)'*sc;  % bin centers
  %for m=1:M, %m = 3; h = histc(A(m,:),b); % explore each channel separately
  h = histc(A(:),b); h = h(:); % all col vecs
  ib = find(h>0.25*max(h)); y = log(h(ib)); x = b(ib); % fit y=const-x^2/(2.eta^2)
  [co] = lscov([1+0*x, x, -x.^2/2],y); eta = 1./sqrt(co(3)); % eta = around 15
  %figure; bar(b,h); set(gca,'yscale','log'); hold on; plot(b,exp(co(1)+co(2)*b-b.^2/(2*eta^2)),'r-');
  noi.eta = eta; % use est std error as Gaussian model

elseif strcmp(meth,'j')      % Jeremy idea: mode of clip-estimated eta distn
  Nc = 1e4;   % # clips
  Tclip = 0.003;  % time clip length in secs. todo: make opt
  T = ceil(Tclip*samplefreq); % # samples in time clip
  [M Nt] = size(A);
  t = randi(Nt-T,1,Nc); % start indices of clips
  etas = nan(1,Nc);
  for c=1:Nc
    etas(c) = sqrt(sum(sum(A(:,t(c):t(c)+T-1).^2)) / (M*T));
  end
  meta = mean(etas); b = meta*(0.5:0.03:10);  % go up to 10x the mean
  f = histc(etas,b);
  %figure; plot(b,f,'+-');
  [~,i] = max(f);      % mode of histogram
  noi.eta = b(i);
end

% todo: meas autocorr & fit to model? This depends on A filtering of course
end