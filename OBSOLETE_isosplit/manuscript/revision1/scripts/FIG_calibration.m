function FIG_calibration

addpath('../../../core');
addpath('../../../core/cluster/isosplit');

close all;

curve_len=200;
fontsize=14;

Ns=[100,1000,10000,100000,1000000];
colors='rgbmc';

figure('Color',[1,1,1]);
set(gcf,'Position',[50,50,1200,400]);
subplot('Position',[0.06,0.13,0.42,0.8]);
for ii=1:length(Ns)
	N=Ns(ii);
	[avg,stdev,scores]=read_isosplit_calibration(N);
	plot(1:curve_len,avg(1:curve_len),colors(ii));
	hold on;
end;
%title('Mean of A_m(\tilde{E}) for various N');
legend('N=100','N=1,000','N=10,000','N=100,000','N=1,000,000','Location','East');
xlabel('m = number of points to average');
ylabel('\mu_{N,m}');
text(0.05,0.92,'A','FontSize',18,'Units','normalize');


subplot('Position',[0.55,0.13,0.42,0.8]);
for ii=1:length(Ns)
	N=Ns(ii);
	[avg,stdev,scores]=read_isosplit_calibration(N);
	plot(1:curve_len,stdev(1:curve_len),colors(ii));
	hold on;
end;
%title('Std. Dev. of A_m(\tilde{E}) for varying N');
legend('N=100','N=1,000','N=10,000','N=100,000','N=1,000,000');
xlabel('m = number of points to average')
ylabel('\sigma_{N,m}');
text(0.05,0.92,'B','FontSize',18,'Units','normalize');

set(gcf,'paperposition',[0,0,8,3]);
print('../calibration_curves_01.eps','-depsc2');

incr=1.02;
ms=[1,5,10,20];
ps=[0.01,0.05,0.1,0.2];
numpts=700;
AVG=zeros(numpts,length(ms));
STDEV=zeros(numpts,length(ms));
SCORE=zeros(numpts,length(ps));
for kk=1:numpts
	N=incr^kk;
	fprintf('.');
	if (mod(kk,25)==0) fprintf('\n'); end;
	[avg,stdev,scores]=read_isosplit_calibration(N);
	for ii=1:length(ms)
		AVG(kk,ii)=avg(ms(ii));
		STDEV(kk,ii)=stdev(ms(ii));
	end;
	for ii=1:length(ps)
		SCORE(kk,ii)=scores(floor(length(scores)*(1-ps(ii))));
	end;
end;
fprintf('\n');

figure('Color',[1,1,1])
set(gcf,'Position',[50,500,1800,400]);
subplot('Position',[0.05,0.13,0.26,0.8]);
for ii=1:length(ms)
	semilogx(incr.^(1:numpts),AVG(:,ii),colors(ii));
	hold on;
end;
xlim([10,incr^(numpts+5)]);
xlabel('N = sample size');
ylabel('\mu_{N,m}');
legend('m = 1','m = 5','m=10','m=20');
text(0.15,0.92,'A','FontSize',18,'Units','normalize');


subplot('Position',[0.37,0.13,0.26,0.8]);
for ii=1:length(ms)
	semilogx(incr.^(1:numpts),STDEV(:,ii),colors(ii));
	hold on;
end;
xlim([10,incr^(numpts+5)]);
xlabel('N = sample size');
ylabel('\sigma_{N,m}');
legend('m = 1','m = 5','m = 10','m = 20');
text(0.1,0.92,'B','FontSize',18,'Units','normalize');

subplot('Position',[0.69,0.13,0.26,0.8]);
for ii=1:length(ps)
	semilogx(incr.^(1:numpts),SCORE(:,ii),colors(ii));
	hold on;
end;
xlim([10,incr^(numpts+5)]);
xlabel('N = sample size');
ylabel('\eta_{N} threshold');
text(0.05,0.92,'C','FontSize',18,'Units','normalize');

legend('\alpha = 0.01','\alpha = 0.05','\alpha = 0.10','\alpha = 0.20');

set(gcf,'paperposition',[0,0,12,3.5]);
print('../calibration_curves_02.eps','-depsc2');

end
