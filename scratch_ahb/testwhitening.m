% basic tests of prewhitening theory. OBSOLETE: see: comparecommonmode.m
% Barnett 5/18/16

if 0 % really basic no-spike tests, iid noise and CM, check analytical model:
Ms = 2.^(2:6);  % range of channel #s
N = 1e5;
alpha = 0.5;     % cm std dev rel to signal std dev
for i=1:numel(Ms), M = Ms(i);
  Y = randn(M,N);               % iid
  cm = randn(1,N);
  Y = Y + repmat(alpha*cm,[M 1]);  % CM, ampl
  [V D] = eig(Y*Y'); D = diag(D); Yw = (V*diag(1./sqrt(D/N))*V')*Y; % jfm
  fprintf('M=%d:\tD/N=[%.3g,%.3g,..] (pred %.3g,%.3g)\tcm rej ratio=%.3g (pred %.3g)\n',M,D(end)/N,mean(D(1:end-1))/N,1+alpha^2*M,1-2*alpha^2,mean(Y*cm')/mean(Yw*cm'),sqrt(1+alpha^2*M))
end
return
% verifies reduction ratio is sqrt(1+alpha^2M).
% So, not much use for small CM signals.
end


addpath ~/validspike/stageA       % generate spikes

clear
M = 16;   % # channels
d.samplerate = 2e4;

K = 5;
Msupport = 3;  % how many elecs the waveform spreads over

T = 30;   % waveform timepts
Tcen = floor((T+1)/2);  % defn of center time of clip
t = (1:T)-Tcen;   % imtes relative to Tcen in timepts
W=zeros(M,T,K);
rng(0);
for k=1:K
  wid = 1+4*rand;
  ampl = 2/wid;
  chans = mod(randi(M)-1 + [0:Msupport-1],M)+1;   % support
  shape = (1.0 + 4.0*(rand-0.5)*t/wid) .* exp(-.5*t.^2/wid^2); % 1st-o poly * gauss
  % same on all supported channels...
  W(chans,:,k) = W(chans,:,k) + ampl*repmat(shape,[Msupport 1]);
end
figure; ms_view_templates(W);

Tsec = 60;  % duration
N = round(Tsec*d.samplerate);         % number of time points
rates = 2.^(0:K-1)/d.samplerate;    % dyadic set of rates
o.amplsig = 0.0;
[times labels ampls] = randomfirings(N,rates,o);
d.truefirings = 'wh_truefirings.mda';
writemda64([0*times;times;labels],d.truefirings);
Y = synthesize_timeseries(W,N,times,labels,ampls);

%Y = 0*Y;   % kill all spikes

cm = 0*0.2*randn(1,N);
o = []; o.freq_min=0; o.freq_max = 6000; o.samplerate = d.samplerate;
cm = ms_bandpass_filter(cm,o);
Y = Y + ones(M,1)*cm;    % outer prod.  same across all chans
eta = 0.2;               % noise std deviation per sample per channel
noi = eta * randn(size(Y));
noi = ms_bandpass_filter(noi,o);
%Y = Y + noi;

Y = Y/sqrt(mean(Y(:).^2));      % make unit variance

d.timeseries = 'wh_raw.mda';
writemda32(Y,d.timeseries);

if 0, thresh = 1.0; o.verb = 0;
  C = empiricalspacetimecorr(struct('A',Y,'dt',1/d.samplerate),thresh,o);
  C(:,:,41)
end % way too slow to plot M^2 bar graphs!!!

%spikespy({Y,times,labels})

%Y*cm'   % how much CM in chans
disp('whiten...'), Yw = ms_whiten(Y);
%Yw = Yw/sqrt(mean(Yw(:).^2));      % make unit variance
%[V D] = eig(Y*Y'); Yw = (diag(1./sqrt(diag(D)/N))*V')*Y;  % 1/N so var stays 1 
% ahb old method, puts it all in last channel

%Yw*cm'   % how much CM in chans

mean(Y(:).^2), mean(Yw(:).^2)  % check rms same

spikespy({Y,'Y'},{Yw,'Yw'});

fprintf('CM proj ratio = %.3g\n',mean(Y*cm')/mean(Yw*cm'))  % 2.9x reduction in projection onto CM
% (why not 4x = sqrt(M) ? )

if 0, [V D] = eig(Y*Y'); %V
  disp('eigvals YYt:')
  diag(D)/N/M  % one eigval 16x the others, with const eigvec
               % (since uncorr noise and CM noise are equal ampl)
  disp('eigvals YwYwt:')
  [V D] = eig(Yw*Yw'); %V
  diag(D)/N/M  % all eigvals equal as expect
end





