% test filter, Matlab vs C.  Barnett 2/16/16

clear; in = [tempdir 'testfilt0.mda']; out = [tempdir 'testfilt1.mda'];
M=8; N=1e7; X = zeros(M,N); X(:,N/2) = 1; %X = randn(M,N); % dummy data spike
writemda(X,in);
samplerate = 2e4;
T = N/samplerate; df = 1/T; f=df*[0:N/2 -N/2+1:-1];
j = 1:round(N/1e4):N/2+1;  % which freq indices to plot
ofilt.samplefreq = samplerate;
ofilt.freq_min=300;
ofilt.freq_max=inf;

mscmd_bandpass_filter(in,out,ofilt); Y = readmda(out);
tic; Ym = ms_filter(X,ofilt); fprintf('ms matlab filter %.3g\n',toc)
fprintf('rel output err btw matlab & C: %.3g\n',norm(Ym-Y)/norm(Y))
Xh = fft(X); Yh = fft(Y);
figure; subplot(2,1,1); t = (1:N)/samplerate; plot(t,X,'+-',t,Y,'o-');
axis([T/2-.002 T/2+.002 -.5 1]); title('impulse response');
subplot(2,1,2); plot(f(j),abs(Xh(:,j)).^2,'-');
hold on; plot(f(j),abs(Yh(:,j)).^2,'r-');
title('unfiltered and filtered |Xhat|'); drawnow
