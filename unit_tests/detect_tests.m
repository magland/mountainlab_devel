function detect_tests(noise_level)

if (nargin<1) noise_level=4; end; %try noise_level=0

close all;

N=400;
num_cycles=40;

%The timeseries
a=linspace(-1,1,N);
X=test_function(a,num_cycles); X=X+randn(size(X))*(noise_level+0.1); %be sure to add at least a bit of noise

%This is just for reference
a2=linspace(-1,1,N*10);
X2=test_function(a,num_cycles)*1.5; 

example_opts.num_cycles=num_cycles;
example_opts.betas = [1,2,10];
opts.detect_threshold=8;
opts.detect_interval=6;
opts.clip_size=10;
opts.sign=1; opts.polarity='p';

example_opts.use_pca=0;
example_opts.use_mscmd=0;
opts.pca_denoise_jiggle=0;
run_detect_test(X,X2,example_opts,opts);

example_opts.use_pca=1;
example_opts.use_mscmd=0;
opts.pca_denoise_jiggle=0;
run_detect_test(X,X2,example_opts,opts);

example_opts.use_pca=0;
example_opts.use_mscmd=1;
opts.pca_denoise_jiggle=0;
run_detect_test(X,X2,example_opts,opts);

example_opts.use_pca=1;
example_opts.use_mscmd=1;
opts.pca_denoise_jiggle=0;
run_detect_test(X,X2,example_opts,opts);

example_opts.use_pca=1;
example_opts.use_mscmd=1;
opts.pca_denoise_jiggle=3;
run_detect_test(X,X2,example_opts,opts);

end

function ret=test_function(a,num_cycles)
ret=cos(a*pi*num_cycles)*20;
end

function run_detect_test(X,X2,example_opts,opts)

N=length(X);
betas=example_opts.betas;
num_cycles=example_opts.num_cycles;
if (example_opts.use_pca)
    opts.num_pca_denoise_components=3;
    opts.meth='p';
else
    opts.num_pca_denoise_components=0;
    opts.meth='x';
end;
    
fA=figure; set(gcf,'position',[100 100 2000 1200]);  % nice wide figure
for i=1:numel(betas)
  if example_opts.use_mscmd
      opts.upsampling_factor = betas(i);
      times=run_mscmd_detect3(X,opts);
  else
      opts.beta = betas(i);
      times=ms_detect3(X,opts);
  end;

  figure(fA);
  subplot(numel(betas),1,i);
  plot(linspace(1,length(X),length(X2)),X2,'g'); hold on;
  plot(1:length(X),X,'k'); hold on;
  for j=1:length(times), plot(times(j)*[1 1],ylim,'r'); end   % times as vlines
  pct=(times-1)/(N-1);
  points_per_cycle=N/num_cycles;
  errs=abs(pct*num_cycles-round(pct*example_opts.num_cycles));
  errs_timepoints=errs*points_per_cycle;
  title0=sprintf('use-mscmd=%d, beta=%d, use-pca=%d, jiggle=%d, err (timepoints): %g',example_opts.use_mscmd,betas(i),example_opts.use_pca,opts.pca_denoise_jiggle,mean(errs_timepoints));
  title(title0);
end;

end

function times=run_mscmd_detect3(X,opts)
path1=[tempdir,'/test1.mda'];
path2=[tempdir,'/test2.mda'];
writemda32(X,path1);
mscmd_detect3(path1,path2,opts);
detect=readmda(path2);
times=detect(2,:);
end