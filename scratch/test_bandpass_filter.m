close all;

M=8;
N=4e7;

%X=zeros(M,N);
%for j=1:size(X,1)
%    X(j,mod(50+j*10,size(X,2))+1)=1;
%end;

X=randn(M,N);

tic;
writemda(X,'scratch/XX.mda');
fprintf('time for writing: %g\n',toc);

cmd='/home/magland/dev/mountainsort/cpp/bin/mountainsort filter --samplefreq=30000 --freq_min=300 --freq_max=10000 --input=/home/magland/dev/mountainsort/scratch/XX.mda --output=/home/magland/dev/mountainsort/scratch/YY.mda';
tic;
system(cmd);
fprintf('time for mountainsort filter: %g\n',toc);

tic;
Y=readmda('scratch/YY.mda');
fprintf('time for reading: %g\n',toc);

tic;
Y2=ms_filter(X,struct('samplefreq',30000,'freq_min',300,'freq_max',10000));
fprintf('time for ms_filter: %g\n',toc);

%figure; plot(1:length(Y),Y,'b',1:length(Y),Y2,'r',1:length(Y),Y-Y2,'k');

Yb=Y(:,500:end-500);
Y2b=Y2(:,500:end-500);

max(abs(Yb(:)-Y2b(:)))
