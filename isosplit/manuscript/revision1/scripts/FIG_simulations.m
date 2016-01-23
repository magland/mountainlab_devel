function FIG_simulations

close all;

addpath('3rdparty/vlfeat-0.9.20/toolbox');
vl_setup; % for gmm using vlfeat

addpath('/home/magland/dev/mountainsort/isosplit');
addpath('internal');
addpath('cluster');

gfopts.num_repeats=20;
gfopts.do_isosplit=1;
gfopts.do_kmeans=1;
gfopts.do_gmm=1;
gfopts.do_gmm_bic=1;
gfopts.do_dbscan=1;
gfopts.show_all_cluster_plots=0;
gfopts.do_3_clusters=1;
gfopts.do_6_clusters=1;
gfopts.do_12_clusters=1;
gfopts.output_basename='simoutput/output';

run_simulation1(gfopts);
run_simulation2(gfopts);
run_simulation3(gfopts)
run_simulation4(gfopts);
run_simulation5(gfopts);
%%%% run_simulation6(gfopts); %noise dimensions

end

function run_simulation1(gfopts)

rng(1); % Use figure 4

opts1.ndims=2;
opts1.zdist=2.5;
opts1.pop_range=[500,500];
opts1.sigma_scale=1;
opts1.spread_factor=0;
opts1.anisotropy_factor=0;
opts1.nongaussian=0;
opts1.verbose=0;

dbeps1=0.27; dbeps2=[];

simopts.num_repeats=gfopts.num_repeats;

if gfopts.do_3_clusters
	opts1.num_clusters=3;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results1=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if gfopts.do_6_clusters
	opts1.num_clusters=6;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results2=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if gfopts.do_12_clusters
	opts1.num_clusters=12;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results3=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

print_results_as_latex_table(algopts,{results1,results2,results3},[gfopts.output_basename,'_sim1.txt']);

end


function run_simulation2(gfopts,nongaussian)

if nargin<2 nongaussian=0; end;

rng(1); % Use figure 4 for simulation2, figure 4 for simulation3

opts1.ndims=2;
opts1.zdist=2.5;
opts1.pop_range=[100,1000];
opts1.sigma_scale=1;
opts1.spread_factor=2;
opts1.anisotropy_factor=1.2;
opts1.nongaussian=nongaussian;
opts1.verbose=0;

dbeps1=0.3; dbeps2=[];

simopts.num_repeats=gfopts.num_repeats;

if gfopts.do_3_clusters
	opts1.num_clusters=3;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results1=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if gfopts.do_6_clusters
	opts1.num_clusters=6;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results2=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if gfopts.do_12_clusters
	opts1.num_clusters=12;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results3=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if ~nongaussian
	print_results_as_latex_table(algopts,{results1,results2,results3},[gfopts.output_basename,'_sim2.txt']);
else
	print_results_as_latex_table(algopts,{results1,results2,results3},[gfopts.output_basename,'_sim3.txt']);
end;

end

function run_simulation3(gfopts)
run_simulation2(gfopts,1);
end

function run_simulation4(gfopts)

rng(1); % Use run 10

opts1.ndims=2;
opts1.zdist=1.7;
opts1.pop_range=[500,500];
opts1.sigma_scale=1;
opts1.spread_factor=0;
opts1.anisotropy_factor=0;
opts1.nongaussian=0;
opts1.verbose=0;

dbeps1=0.3; dbeps2=[];

simopts.num_repeats=gfopts.num_repeats;

if gfopts.do_3_clusters
	opts1.num_clusters=3;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results1=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if gfopts.do_6_clusters
	opts1.num_clusters=6;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results2=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if gfopts.do_12_clusters
	opts1.num_clusters=12;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results3=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

print_results_as_latex_table(algopts,{results1,results2,results3},[gfopts.output_basename,'_sim4.txt']);

end

function run_simulation5(gfopts)

rng(1); % Use figure 19

opts1.ndims=6;
opts1.zdist=2.5;
opts1.pop_range=[100,1000];
opts1.sigma_scale=1;
opts1.spread_factor=2;
opts1.anisotropy_factor=1.2;
opts1.nongaussian=0;
opts1.verbose=0;

dbeps1=1.5; dbeps2=[];

simopts.num_repeats=gfopts.num_repeats;

if gfopts.do_3_clusters
	opts1.num_clusters=3;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results1=run_simulation(gfopts,opts1,algs,algopts,simopts);
else
	results1=[];
end;

if gfopts.do_6_clusters
	opts1.num_clusters=6;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results2=run_simulation(gfopts,opts1,algs,algopts,simopts);
else
	results2=[];
end;

if gfopts.do_12_clusters
	opts1.num_clusters=12;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results3=run_simulation(gfopts,opts1,algs,algopts,simopts);
else
	results3=[];
end;

print_results_as_latex_table(algopts,{results1,results2,results3},[gfopts.output_basename,'_sim5.txt']);

end

function run_simulation6(gfopts)

rng(1); % Use figure 4 for simulation2, figure 4 for simulation3

opts1.ndims=2;
opts1.zdist=2.5;
opts1.pop_range=[100,1000];
opts1.sigma_scale=1;
opts1.spread_factor=2;
opts1.anisotropy_factor=1.2;
opts1.nongaussian=0;
opts1.verbose=0;

dbeps1=0.3; dbeps2=[];

simopts.num_repeats=gfopts.num_repeats;
simopts.num_noise_dimensions=30;

if gfopts.do_3_clusters
	opts1.num_clusters=3;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results1=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if gfopts.do_6_clusters
	opts1.num_clusters=6;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results2=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if gfopts.do_12_clusters
	opts1.num_clusters=12;
	[algs,algopts]=get_algorithms(gfopts,opts1.num_clusters,dbeps1,dbeps2);
	results3=run_simulation(gfopts,opts1,algs,algopts,simopts);
end;

if ~nongaussian
	print_results_as_latex_table(algopts,{results1,results2,results3},[gfopts.output_basename,'_sim2.txt']);
else
	print_results_as_latex_table(algopts,{results1,results2,results3},[gfopts.output_basename,'_sim3.txt']);
end;

end



function results=run_simulation(gfopts,opts1,algs,algopts,simopts)

if (~isfield(simopts,'num_noise_dimensions')) simopts.num_noise_dimensions=0; end;

nr=2; nc=3;
fontsize=14;

subplot_opts.ymargin=0.04;
subplot_opts.yspacing=0.06;
subplot_opts.xmargin=0.03;
subplot_opts.xspacing=0.03;

accuracies=zeros(length(algs),simopts.num_repeats);
times=zeros(length(algs),simopts.num_repeats);

for rr=1:simopts.num_repeats

if (~isfield(opts1,'pops'))
	opts1.pops=randi(opts1.pop_range,1,opts1.num_clusters);
end;

fprintf('########## rr=%d\n',rr);	

[clusters,samples,labels]=generate_random_clusters(opts1);
	
%if (rr==1)||(gfopts.show_all_cluster_plots)
%	HFF=figure; set(gcf,'Position',[50,50,1200,800]);
%	subplot_jfm(nr,nc,1,subplot_opts);
	%plot_samples(samples,labels,10); title('True');
%end;

for j=1:length(algs)
	alg=algs{j};
	opts=algopts{j};
	
	if (simopts.num_noise_dimensions>0)
		samples=cat(1,samples,1*randn(simopts.num_noise_dimensions,size(samples,2)));
	end;
	
	tA=tic;
	labels0=alg(samples,opts);
	elapsed=toc(tA);
	times(j,rr)=elapsed;
	
	if (simopts.num_noise_dimensions>0)
		samples=samples(1:opts1.ndims,:);
	end;
	
	[~,labels0]=ss_getbestshuffling(labels0,labels);
	
	acc0=mean(compute_accuracies(labels,labels0));
	str=opts.name; while (length(str)<30) str=[str,' ']; end;
	fprintf('%s\tAccuracy = %.1f\tElapsed = %g\n',str,acc0*100,elapsed);
	accuracies(j,rr)=acc0;
	
% 	if (rr==1)||(gfopts.show_all_cluster_plots)
% 		figure(HFF);
% 		subplot_jfm(nr,nc,j+1,subplot_opts);
% 		str=opts.name;
% 		str=strrep(str,'$','');
% 		plot_samples(samples,labels0,10); title(str);
% 		drawnow;
% 	end;
end;
fprintf('\n');

end;

fprintf('\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
results=zeros(length(algs),3);
for j=1:length(algs)
	alg=algs{j};
	opts=algopts{j};
	str=opts.name; while (length(str)<30) str=[str,' ']; end;
	avg0=mean(accuracies(j,:))*100;
	std0=sqrt(var(accuracies(j,:)*100))/sqrt(simopts.num_repeats);
	time0=mean(times(j,:));
	fprintf('%s\tavg.acc=%.1f\tstd.acc=%.2f\tavg.time=%g\n',str,avg0,std0,time0);
	results(j,1)=avg0;
	results(j,2)=std0;
	results(j,3)=time0;
end;
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');

end

function print_results_as_latex_table(algopts,results,fname)

str='';
for j=1:length(algopts)
	str=[str,sprintf('\t%s',algopts{j}.name)];
	for k=1:length(results)
		result0=results{k};
		if (length(result0)>0)
			str=[str,sprintf(' & $%.1f \\pm %.1f\\%%$',result0(j,1),result0(j,2))];
		else
			str=[str,sprintf(' & $0 \\pm 0\\%%$')];
		end;
	end;
	str=[str,sprintf(' \\\\\n')];
end;

str=[str,sprintf('\n')];
for j=1:length(algopts)
	str=[str,sprintf('\t%s',algopts{j}.name)];
	for k=1:length(results)
		result0=results{k};
		if (length(result0)>0)
			str=[str,sprintf(' & $%.2f$',result0(j,3))];
		else
			str=[str,sprintf(' & $0$')];
		end;
	end;
	str=[str,sprintf(' \\\\\n')];
end;

fprintf('\n\n%s',str);
FF=fopen(fname,'w');
fprintf(FF,'%s',str);
fclose(FF);

end

function simulation_result_graph

str=fileread('../iso-split.tex');
lines=strsplit(str,'\n');
DATA1=zeros(0,3);
DATA2=zeros(0,3);
for j=1:length(lines)
	line=lines{j};
	if (length(strfind(line,'\pm'))>0)&&(length(strfind(line,' & '))>2)
		disp(line);
		line=strrep(line,'$','');
		line=strrep(line,' ','');
		line=strrep(line,'\%','');
		line=strrep(line,'\\','');
		strings=strsplit(line,'&');
		tmp1=[];
		tmp2=[];
		for k=2:length(strings)
			str=strings{k};
			vals=strsplit(str,'\\pm');
			tmp1(end+1)=eval(vals{1});
			tmp2(end+1)=eval(vals{2});
		end;
		DATA1(end+1,:)=tmp1;
		DATA2(end+1,:)=tmp2;
	end;
end;

DATA1=reshape(DATA1,5,size(DATA1,1)/5,size(DATA1,2));
DATA2=reshape(DATA2,5,size(DATA2,1)/5,size(DATA2,2));

disp(DATA1);
disp(DATA2);

for cc=1:3
	figure('name',sprintf('cc=%d',cc)); set(gcf,'Color','w');
	barweb(DATA1(:,:,cc)',DATA2(:,:,cc)'); hold on;
	xlabel('Simulation');
	ylabel('Avg. Accuracy (%)');
	legend('ISO-SPLIT','K-means (K known)','GMM (K known)','GMM (BIC)','DBSCAN');
end;

end

function plot_samples(samples,labels,marker_size)

if (nargin<3) marker_size=1; end;

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

function [acc,false_split,false_merge,unclassified]=compute_accuracies(labels1,labels2)

K1=max(labels1);
K2=max(labels2);
acc=zeros(K1,1);
false_split=zeros(K1,1);
false_merge=zeros(K1,1);
for k1=1:K1
	best_choice=1;
	best_count=0;
	for k2=1:K2
		ct=length(find((labels1==k1)&(labels2==k2)));
		if (ct>best_count)
			best_choice=k2;
			best_count=ct;
		end;
	end;
	k2=best_choice;
	numer=length(find((labels2>0)&((labels1==k1)&(labels2==k2))));
	denom1=length(find((labels2>0)&(labels1==k1)));
	denom2=length(find((labels2>0)&(labels2==k2)));
	if (denom1==0) denom1=1; end;
	if (denom2==0) denom2=1; end;
	acc(k1)=min(numer/denom1,numer/denom2);
	
	numer=length(find((labels2>0)&((labels1==k1)&(labels2==k2))));
	denom=length(find((labels2>0)&((labels1==k1))));
	if (denom>0)
		false_split(k1)=(1-numer/denom);
	end;
	
	numer=length(find((labels2>0)&((labels1==k1)&(labels2==k2))));
	denom=length(find((labels2>0)&((labels2==k2))));
	if (denom>0)
		false_merge(k1)=(1-numer/denom);
	end;
end;
unclassified=length(find(labels2<=0))/length(labels2);

end
