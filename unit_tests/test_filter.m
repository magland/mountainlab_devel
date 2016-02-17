% test filter, Matlab vs C.  Barnett 2/16/16

clear; in = [tempdir 'testfilt0.mda']; out = [tempdir 'testfilt1.mda'];
M=8; N=1e7;
X = zeros(M,N); X(:,N/2) = 1;    % dummy data spike
%X = randn(M,N);                 % white noise
writemda(X,in);
samplerate = 2e4;
T = N/samplerate; df = 1/T; f=df*[0:N/2 -N/2+1:-1];
ofilt.samplefreq = samplerate;
ofilt.freq_min=300;
ofilt.freq_max=inf;

mscmd_bandpass_filter(in,out,ofilt); Y = readmda(out);
tic; Ym = ms_filter(X,ofilt); fprintf('ms matlab filter %.3g\n',toc)
fprintf('rel output err btw matlab & C: %.3g\n',norm(Ym-Y)/norm(Y))
figure; subplot(2,1,1); t = (1:N)/samplerate;
trng = T/2 + [-1 1]*.002; i = (t>trng(1) & t<trng(2));  % t indices to plot
plot(t(i),X(:,i),'+-',t(i),Y(:,i),'o-');
axis([trng -.5 1]); title('impulse response');
j = 1:round(N/1e4):N/2+1;                               % freq indices to plot
Xh = fft(X(1,:)); Yh = fft(Y(1,:));                     % just check first row
subplot(2,1,2); plot(f(j),abs(Xh(j)).^2,'-');
hold on; plot(f(j),abs(Yh(j)).^2,'r-');
title('unfiltered and filtered first row of |Xhat|'); drawnow
