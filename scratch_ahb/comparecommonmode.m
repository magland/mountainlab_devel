function comparecommonmode
% Compare common mode subtraction methods, via ratio of proj onto CM after/before
% Needs mountainlab in path, and spikespy.
% Barnett 5/23/16

d.samplerate = 2e4;
M = 16;
Tsec = 60;    % duration
N = round(Tsec*d.samplerate);         % number of time points
K = 5; Msupport = 3; T = 30;  % if spikes, how many, channel support, clip wid
rng(0);          % repeatable expts
W = synthwaveforms(M,Msupport,K,T); %figure; ms_view_templates(W); drawnow
rates = 10.0 * ones(1,K)/d.samplerate;   % all 10 Hz
%2.^(0:K-1)/d.samplerate;    % dyadic set of rates, not dominating vs noi
Asig = 1.0;          % ampl height of spikes rel to iid noise only (not CM)
Anoi = 0.2;  % 0.05
o.amplsig = 0.0;
[times labels ampls] = randomfirings(N,rates,o);  % fix once and for all
fprintf('synth %d spikes with K=%d\n',numel(times),K)
ampls = ampls * Asig;     % how big to make signal
Ysig = synthesize_timeseries(W,N,times,labels,ampls);
fprintf('mean sq of sig: %.3g\n',mean(Ysig(:).^2))
Ynoi = randn(M,N);               % iid noise
cm = randn(1,N);                 % common mode iid "noise"
if 1                    % whether to bandpass filter the noise
  o = []; o.freq_min=0; o.freq_max = 6000; o.samplerate = d.samplerate;  % filter
  Ynoi = ms_bandpass_filter(Ynoi,o); cm = ms_bandpass_filter(cm,o);
  Ynoi = Ynoi*(1/sqrt(mean(Ynoi(:).^2))); cm = cm*(1/sqrt(mean(cm.^2))); % var=1
end
Ynoi = Anoi * Ynoi;
fprintf('mean sq of noi: %.3g\n',mean(Ynoi(:).^2))

meth{1} = 'mean subtract'; f{1} = @(Y) subcm(Y,mean(Y,1));
meth{2} = 'jfm singvals->1'; f{2} = @(Y) ms_whiten(Y);  % names and methods
meth{3} = 'axel trimmin'; f{3} = @(Y) subcm(Y,trimmin(Y,60,M));  % 1-60%=6/16
meth{4} = 'ahb wei trimmin'; f{4} = @(Y) subcm(Y,estcm(Y,0.4));  % override lev
meths = 1:4;   % which to test

for spikes = 0; %[0 1]       % ================= two expts
  als = 0:0.05:0.4;     % cm std dev rel to other iid noise std dev
  cmrs = nan(max(meths),numel(als)); prs=cmrs; % cm rej ratios, pred ratio, <=1
  for i=1:numel(als), al = als(i)   % ---------- loop over cm noise levels
    fprintf('mean sq of CM: %.3g\n',al^2*mean(cm.^2))
    Y = spikes*Ysig + Ynoi + repmat(al*cm,[M 1]);     % forward model
    [V D] = eig(Y*Y'); D = diag(D);
    fprintf('max lambda/N = %.3g (pred %.3g),\tmean of other lambda/N = %.3g (pred %.3g)\n',D(end)/N,Anoi^2+al^2*M,mean(D(1:end-1))/N,Anoi^2)  % good pred when no spikes
    Y = Y*(1/sqrt(mean(Y(:).^2)));   % make unit variance overall
    cmip = norm(Y*cm');   % common mode inner prod in raw data (Y not var=1)
    for m=meths    % do methods (on var=1 data) and metrics
      Yw = f{m}(Y);
      Yw = Yw*(1/sqrt(mean(Yw(:).^2)));   % make unit variance overall, fairness
      cmrs(m,i) = norm(Yw*cm') / cmip;
      if i==numel(als) && spikes, Ywkeep{m} = Yw; end % keep the max alpha cases
    end
    % some predictions
    prs(2,i) = Anoi/sqrt(Anoi^2+al^2*M);
    [cmrs(meths,i)',prs(meths,i)']  % show current metrics
  end                               % -----------
  % spikes=0 case verifies reduction is sqrt(1+(alpha/Anoi)^2.M)
  figure; plot(als,cmrs(meths,:),'+-'); hold on;
  plot(als,prs(2,:),'o-');   % pred only for meth 2
  legs = {meth{meths}, [meth{2}, ' pred']}; legend(legs);
  xlabel('\alpha (CM ampl)'); ylabel('\rho (CM inner-prod ratio)');
  title(sprintf('spikes=%d, Anoi=%g: CM inner-prod ratio \\rho vs \\alpha',spikes,Anoi))
  drawnow
end
%print -depsc2 ~/ss/comparecmrho_Anoi02.eps

if 0
% compare meth outputs w/spikes at maximal alpha case...
addpath ~/spikespy/matlab/  % use old spikespy
spikespy({Y,times,labels,'Y'},{Ywkeep{1},['Yw ',meth{1}]},{Ywkeep{2},['Yw ',meth{2}]},{Ywkeep{3},['Yw ',meth{3}]},{Ywkeep{4},['Yw ',meth{4}]});

figure; subplot(1,numel(meths)+1,1);
Wo = ms_templates(ms_extract_clips2(Y,times,T),labels); ms_view_templates(Wo);
fprintf('raw: mean template peak range = [%.3g,%.3g]\n',min(squeeze(max(max(abs(Wo),[],1),[],2))),max(Wo(:)))
for m=meths, subplot(1,numel(meths)+1,m+1);  % show mean W's using true times...
  Ww = ms_templates(ms_extract_clips2(Ywkeep{m},times,T),labels);
  fprintf('meth=%d: mean template peak range = [%.3g,%.3g]\n',m,min(squeeze(max(max(abs(Ww),[],1),[],2))),max(Ww(:)))
  ms_view_templates(Ww);
  title(['mean W, ',meth{m}])
end
set(gcf,'paperposition',[0 0 15 10]); print -depsc2 meanW.eps
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Yw = subcm(Y,cm)           % subtract a common mode sig from all chan
M = size(Y,1);
Yw = Y - repmat(cm,[M 1]);

function cm = estcm(Y,lev)           % ahb weighted mean, nonlin like trimmin
% but not discont like trimmin is. lev (opt) = amplitude scale where wei dies.
if nargin<2, lev = 1.0 * sqrt(mean(Y(:).^2)); end
[M N] = size(Y);
wei = exp((-.5/lev^2)*Y.^2);      % Gaussian-decaying weights per signal point
cm = sum(wei.*Y,1) ./ sum(wei,1);  % normalized by col-sums
fprintf('\t estcm: mean sum of weights = %.3g\n',mean(sum(wei,1)))

function trimminTrace = trimmin(data,PCTexclude,numSites,TMchunkSize)
% used by Carl Schoonover at Axel lab, May 2016

if nargin<4, TMchunkSize = 1e7; end                      % ahb added
numToInclude = round((100-PCTexclude)/100 * numSites);   % they said about 6
trimminTrace = zeros(1,length(data));
if TMchunkSize < length(data)
  numChunks = floor(length(data)/TMchunkSize) + 1;
else
  numChunks = 1;
end
% Loop approach can do all-in-one or chunk (for RAM-limited platforms)
for j = 1:numChunks
  chunkBins = (j-1)*TMchunkSize+1 : min([j*TMchunkSize length(data)]);
  dataChunk = data(:,chunkBins);
  [~,I] = sort(abs(dataChunk),'ascend');
  keepIndex = I(1:numToInclude,:) + ones(numToInclude,1) * numSites * (0:length(dataChunk)-1);
  trimminTrace(chunkBins) = mean(dataChunk(keepIndex));
  %disp(['Trimmin chunk ' int2str(j) ' of ' int2str(numChunks) '.']);
end

function W = synthwaveforms(M,Msupport,K,T)
% W = synthwaveforms(M,Msupport,K,T)
%
% Barnett 5/23/16
Tcen = floor((T+1)/2);  % defn of center time of clip
t = (1:T)-Tcen;   % imtes relative to Tcen in timepts
W=zeros(M,T,K);
for k=1:K
  wid = 1+4*rand;
  chans = mod(randi(M)-1 + [0:Msupport-1],M)+1;   % support
  shape = (1.0 + 4.0*(rand-0.5)*t/wid) .* exp(-.5*t.^2/wid^2); % 1st-o poly * gauss
  ampl = 1/max(shape);   % unit peak ampl
  %ampl = 2/wid;
  % same on all supported channels...
  W(chans,:,k) = W(chans,:,k) + ampl*repmat(shape,[Msupport 1]);
end
%disp('waveform peaks before SNR-scaling:'), squeeze(max(max(abs(W),[],1),[],2))'
