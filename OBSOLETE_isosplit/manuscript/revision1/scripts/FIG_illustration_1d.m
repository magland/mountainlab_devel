function FIG_illustration_1d

addpath('../../../core');
addpath('../../../core/cluster/isosplit');

close all;

rng(1);

N=1000;

fontsize=16; % For A,B,C,...

figure;
set(gcf,'Color','w');
set(gcf,'Position',[100,100,1400,900]);

letters={};
letters{1}.a='A'; letters{1}.b='D'; letters{1}.c='G';
letters{2}.a='B'; letters{2}.b='E'; letters{2}.c='H';
letters{3}.a='C'; letters{3}.b='F'; letters{3}.c='I';

subplot_opts.xmargin=0.04;
subplot_opts.ymargin=0.04;
subplot_opts.xspacing=0.04;
subplot_opts.yspacing=0.04;

for pass=1:3
	if (pass==1)
		X=randn(1,N)*0.9;
	elseif (pass==2)
		X=randn(1,N/2)*0.6-1.4;
		X=[X,randn(1,N/2)*0.6+1.4];
	elseif (pass==3)
		X=(rand(1,N)-0.5)*5.5;
		X=X(abs(X)>0.15);
	end;
	xx=[-3,3];
	yy=[-4,4];

	X=sort(X);

	subplot_jfm(3,3,1+(pass-1),subplot_opts);
	hist(X,200);
	xlim(xx);
	set(gca,'xtick',[]);
	if (pass==1) ylabel('Number of samples'); end;
	text(0.05,0.9,letters{pass}.a,'FontSize',fontsize,'Units','normal');
	ylim([0,20]);

	spacings=X(2:end)-X(1:end-1);
	XX=(X(1:end-1)+X(2:end))/2;
	densities=1./spacings;
	D=log(densities);

	Dfit=jisotonic(D,'updown');
	subplot_jfm(3,3,4+(pass-1),subplot_opts);
	h=plot(XX,D,'b',XX,Dfit,'k');
	set(h(2),'linewidth',3);
	xlim(xx);
	set(gca,'xtick',[]);
	if (pass==1) ylabel('Log density'); end;
	text(0.05,0.9,letters{pass}.b,'FontSize',fontsize,'Units','normal');

	Dresid=D-Dfit;
	Dresidfit=jisotonic(Dresid,'downup');
	subplot_jfm(3,3,7+(pass-1),subplot_opts);
	h=plot(XX,Dresid,'b',XX,Dresidfit,'k');
	set(h(2),'linewidth',3);
	xlim(xx); ylim(yy);
	set(gca,'xtick',[]);
	if (pass==1) ylabel('Log density residual'); end;
	text(0.05,0.9,letters{pass}.c,'FontSize',fontsize,'Units','normal');
end;

set(gcf,'paperposition',[0,0,12,7]);
print('../illustration_1d.eps','-depsc2');

end

