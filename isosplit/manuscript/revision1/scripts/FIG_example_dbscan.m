function FIG_example_dbscan

addpath('../../../core');

close all;
rng(10);
list={};

dx0=12;
dy0=8;
dbeps1=0.15;
dbeps2=0.30;

N0=4000;
theta=rand(1,N0)*pi;
x0=10*cos(theta); y0=10*sin(theta);
A0=cat(1,x0,y0); A0=A0+randn(size(A0))*1;
list{end+1}=A0;

N0=2000;
theta=rand(1,N0)*pi;
x0=cat(2,zeros(1,N0/2)-0.5,zeros(1,N0/2)+2); y0=cat(2,zeros(1,N0/2)+3,zeros(1,N0/2)+3);
A0=cat(1,x0,y0); A0=A0+randn(size(A0))*1;
list{end+1}=A0;

A=zeros(2,0);
labels=zeros(1,0);
for j=1:length(list)
	A=cat(2,A,list{j});
	labels=cat(2,labels,ones(1,size(list{j},2))*j);
end;

true_K=max(labels);

gfopts.do_isosplit=1;
gfopts.do_kmeans=1;
gfopts.do_gmm=1;
gfopts.do_gmm_bic=0;
gfopts.do_dbscan=1;
[algs,algopts]=get_algorithms(gfopts,true_K,dbeps1,dbeps2);

run_example(A,labels,algs,algopts);

set(gcf,'paperposition',[0,0,9,6]);
print('../example_dbscan.eps','-depsc2');

end

function run_example(A,labels,algs,algopts)

N=size(A,2);
N

nr=2; nc=3;
fontsize=14;
HFF=figure; set(gcf,'Position',[50,50,1200,800]);

subplot_opts.ymargin=0.04;
subplot_opts.yspacing=0.06;
subplot_opts.xmargin=0.03;
subplot_opts.xspacing=0.03;

figure(HFF);
subplot_jfm(nr,nc,1,subplot_opts);
plot_samples(A,labels,0.1); title('True');

for j=1:length(algs)
	alg=algs{j};
	opts=algopts{j};
	
	disp(opts.name);
	tA=tic;
	labels0=alg(A,opts);
	elapsed=toc(tA);
	fprintf('Elapsed = %g\n',elapsed);
	
	[~,labels0]=ss_getbestshuffling(labels0,labels);
	
	%fprintf('Accuracy = %g\n',mean(compute_accuracies(labels,labels0)));
	
	figure(HFF);
	subplot_jfm(nr,nc,j+1,subplot_opts);
	str=opts.name;
	str=strrep(str,'$','');
	plot_samples(A,labels0,0.1); title(str);
end;

end

function plot_samples(samples,labels,marker_size)

if (nargin<3) marker_size=4; end;

CC=distinguishable_colors(max(labels),{'w'});
colors={};
for j=1:size(CC,1)
	colors{j}=CC(j,:);
end;

for j=1:max(labels)
	xx=samples(1,find(labels==j));
	yy=samples(2,find(labels==j));
	col=colors{mod(j-1,length(colors))+1};
	plot(xx,yy,'.','Color',col,'MarkerSize',marker_size);
	if (j==1) hold on; end;
end;
set(gca,'xtick',[],'ytick',[])

xmin=min(samples(1,:)); xmax=max(samples(1,:));
ymin=min(samples(2,:)); ymax=max(samples(2,:));
xlim([floor(xmin-1),ceil(xmax+1)]);
ylim([floor(ymin-1),ceil(ymax+1)]);

XL=xlim; YL=ylim;
dx=XL(2)-XL(1);
dy=YL(2)-YL(1);
dd=max(dx*1.05,dy*1.05);

xlim([XL(1)-(dd-dx)/2,XL(2)+(dd-dx)/2]);
ylim([YL(1)-(dd-dy)/2,YL(2)+(dd-dy)/2]);

end
