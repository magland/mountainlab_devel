function FIG_decision_boundaries

addpath('../../../core');

close all;

rng(1);

fontsize=16;

simoptions.centers={[0,0],[3.5,4.2],[-1,4.8]};
simoptions.pops={1000,400,400};
simoptions.shapes={[2,2,0],[0.5,0.5,0],[1.5,0.4,0]};

[X,labels]=generate_samples_decision_boundaries(simoptions);

labels1=ss_cluster(X,struct('cmethod','iso'));
labels2=ss_cluster(X,struct('cmethod','k++','K',3));

figure;
set(gcf,'Position',[50,50,1600,500]);
subplot_opts.yoffset=-0.02;

subplot_jfm(1,3,1,subplot_opts);
plot_samples_2d(X,labels);
title('Truth');

subplot_jfm(1,3,2,subplot_opts);
plot_samples_2d(X,labels1); hold on;
pt1=[1.9,4.0];
pt2=[-5.2,2.9]; pt2=pt2+(pt2-pt1)*10;
pt3=[2.4,5.1]; pt3=pt3+(pt3-pt1)*10;
pt1b=[1.7,4.1];
pt4=[4.2,1.2]; pt4=pt4+(pt4-pt1b)*10;
plot([pt1b(1),pt2(1)],[pt1b(2),pt2(2)],'k--');
plot([pt1(1),pt3(1)],[pt1(2),pt3(2)],'k--');
plot([pt1b(1),pt4(1)],[pt1b(2),pt4(2)],'k--');
title('ISO-SPLIT');

subplot_jfm(1,3,3,subplot_opts);
plot_samples_2d(X,labels2); hold on; 
pt1=[0.8,2.1];
pt2=[-5.5,0.5]; pt2=pt2+(pt2-pt1)*10;
pt3=[1.0,5.7]; pt3=pt3+(pt3-pt1)*10;
pt4=[5.5,-1.2]; pt4=pt4+(pt4-pt1)*10;
plot([pt1(1),pt2(1)],[pt1(2),pt2(2)],'k--');
plot([pt1(1),pt3(1)],[pt1(2),pt3(2)],'k--');
plot([pt1(1),pt4(1)],[pt1(2),pt4(2)],'k--');
title('K-means');

set(gcf,'paperposition',[0,0,8,3]);
print('../decision_boundaries.eps','-depsc2');

end

function [samples,labels]=generate_samples_decision_boundaries(opts)

centers=opts.centers;
pops=opts.pops;
shapes=opts.shapes;

samples=zeros(2,0);
labels=zeros(1,0);

for j=1:length(centers)
	xx=randn(1,pops{j});
	yy=randn(1,pops{j});
	shape=shapes{j};
	xx2=xx*shape(1)+yy*shape(3);
	yy2=-xx*shape(3)+yy*shape(2);
	center=centers{j};
	xx2=xx2+center(1);
	yy2=yy2+center(2);
	tmp=zeros(2,pops{j});
	tmp(1,:)=xx2; tmp(2,:)=yy2;
	samples=[samples,tmp];
	labels=[labels,ones(1,pops{j})*j];
end;

end

function plot_samples_2d(samples,labels)

if (nargin<4) subtitle=''; end;

colors='rgbkymc';
for j=1:max(labels)
	xx=samples(1,find(labels==j));
	yy=samples(2,find(labels==j));
	col=colors(mod(j-1,length(colors))+1);
	plot(xx,yy,['.',col]);
	if (j==1) hold on; end;
end;
set(gca,'xtick',[],'ytick',[])

xmin=min(samples(1,:)); xmax=max(samples(1,:));
ymin=min(samples(2,:)); ymax=max(samples(2,:));
xlim([floor(xmin-2),ceil(xmax+2)]);
ylim([floor(ymin-2),ceil(ymax+2)]);

end
