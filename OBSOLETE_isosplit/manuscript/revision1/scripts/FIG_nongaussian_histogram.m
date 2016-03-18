function FIG_nongaussian_histogram

addpath('../../../core');

Z=randn(1,5e7);
Z2=log(abs(Z+3));
figure; hist(Z2,10000);
xlim([-2.5,2.5]);
set(gcf,'Color','w');
set(gca,'ytick',[])

set(gcf,'paperposition',[0,0,5,3.5]);
print('../nongaussian_histogram.eps','-depsc2');

end
