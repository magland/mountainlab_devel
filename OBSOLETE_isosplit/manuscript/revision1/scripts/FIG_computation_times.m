function FIG_computation_times

addpath('../../../core');
addpath('../../../core/cluster/isosplit');

close all;

simopts.num_repeats=50;

 %simulation_time_increasing_N(simopts);
 %simulation_time_increasing_num_clusters(simopts);
 %simulation_time_increasing_ndims(simopts);

%simulation_time_create_figure;

%simulation_alpha_dependence(simopts);
 simulation_K_dependence(simopts);
% simulation_m_max_dependence(simopts);

end

function simulation_time_increasing_N(simopts)

Ns=250:250:5000;
num_repeats=simopts.num_repeats;
num_clusters=6;
K=24;

times=zeros(num_repeats,length(Ns));

if (~exist('simulation_time_increasing_N.mat'))

figure;
for jj=1:length(Ns)
for rr=1:num_repeats
	N=Ns(jj);
	
	fprintf('%d: ',N);
	
	opts.ndims=2;
	opts.zdist=2.5;
	NN=floor(N/6);
	opts.pops=[NN,NN,NN,NN,NN,NN];
	opts.sigma_scale=1;
	opts.num_clusters=num_clusters;
	opts.spread_factor=2;
	opts.anisotropy_factor=1.2;
	opts.num_intersection_trials=1000;
	
	[clusters,samples,labels]=generate_random_clusters(opts);
	tA=tic;
	isosplit(samples,struct('K',24));
	elapsed=toc(tA);
	fprintf('%.4f s\n',elapsed);
	times(rr,jj)=elapsed;
	plot(Ns,mean(times,1),'b.-'); drawnow;
end;
end;

save('simulation_time_increasing_N.mat','Ns','times');

end;

load('simulation_time_increasing_N.mat');

coeffs=polyfit(Ns,mean(times,1),1);
fprintf('coeffs = %.8f, %.8f\n',coeffs(1),coeffs(2));
times_fit=polyval(coeffs,Ns);
figure; plot(Ns,mean(times,1),'b.-',Ns,times_fit,'r');

end

function simulation_time_increasing_num_clusters(simopts)

vals=1:1:20;
N=500;
K=60;
num_repeats=simopts.num_repeats;

if (~exist('simulation_time_increasing_num_clusters.mat'))

times=zeros(num_repeats,length(vals));

figure;
for jj=1:length(vals)
for rr=1:num_repeats
	num_clusters=vals(jj);
	
	fprintf('%d, %d: ',num_clusters,K);
	
	opts.ndims=2;
	opts.zdist=2.5;
	NN=floor(N/num_clusters);
	opts.pops=ones(1,num_clusters)*NN;
	opts.sigma_scale=1;
	opts.num_clusters=num_clusters;
	opts.spread_factor=2;
	opts.anisotropy_factor=1.2;
	opts.num_intersection_trials=1000;
	
	[clusters,samples,labels]=generate_random_clusters(opts);
	tA=tic;
	isosplit(samples,struct('K',K));
	elapsed=toc(tA);
	fprintf('%.4f s\n',elapsed);
	times(rr,jj)=elapsed;
	
	plot(vals,mean(times,1),'b.-'); drawnow;
end;
end;

save('simulation_time_increasing_num_clusters.mat','vals','times');

end;

load('simulation_time_increasing_num_clusters.mat');

coeffs=polyfit(vals,mean(times,1),1);
%fprintf('coeffs = %.8f, %.8f, %.8f\n',coeffs(1),coeffs(2),coeffs(3));
fprintf('coeffs = %.8f, %.8f\n',coeffs(1),coeffs(2));
times_fit=polyval(coeffs,vals);
figure; plot(vals,mean(times,1),'b.-',vals,times_fit,'r');

end

function simulation_time_increasing_ndims(simopts)

vals=1:1:15;
N=500;
K=24;
num_clusters=6;
num_repeats=simopts.num_repeats;

if (~exist('simulation_time_increasing_ndims.mat'))

times=zeros(num_repeats,length(vals));
iteration_counts=zeros(num_repeats,length(vals));

F1=figure; title('Time');
F2=figure; title('Num. Dimensions');
for jj=1:length(vals)
for rr=1:num_repeats
	num_clusters=vals(jj);
	
	fprintf('%d, %d: ',num_clusters,K);
	
	%opts.ndims=vals(jj);
	opts.ndims=2;
	opts.zdist=2.5;
	NN=floor(N/num_clusters);
	opts.pops=ones(1,num_clusters)*NN;
	opts.sigma_scale=1;
	opts.num_clusters=num_clusters;
	opts.spread_factor=2;
	opts.anisotropy_factor=1.2;
	opts.num_intersection_trials=1000;
	
	[clusters,samples,labels]=generate_random_clusters(opts);
	samples=cat(1,samples,randn(opts.ndims-2,size(samples,2)));
	tA=tic;
	[~,info]=isosplit(samples,struct('K',K));
	elapsed=toc(tA);
	fprintf('%.4f s\n',elapsed);
	times(rr,jj)=elapsed;
	iteration_counts(rr,jj)=info.num_iterations;
	
	figure(F1);
	plot(vals,mean(times,1),'b.-'); drawnow;
	figure(F2);
	plot(vals,mean(iteration_counts,1),'r.-'); drawnow;
end;
end;

save('simulation_time_increasing_ndims.mat','vals','times','iteration_counts');

end;

load('simulation_time_increasing_ndims.mat');

coeffs=polyfit(vals,mean(times,1),1);
fprintf('coeffs = %.8f, %.8f\n',coeffs(1),coeffs(2));
times_fit=polyval(coeffs,vals);
figure; plot(vals,mean(times,1),'b.-',vals,times_fit,'r');

end

function simulation_time_create_figure

tmp1=load('simulation_time_increasing_N.mat');
tmp2=load('simulation_time_increasing_num_clusters.mat');
tmp3=load('simulation_time_increasing_ndims.mat');

figure; rr=1; cc=3;
set(gcf,'Position',[50,50,2000,600]);
set(gcf,'Color','w');
ystr='Computation time (ms)';
set(gcf,'Units','normal');
font_size=16;

subplot_opts.xspacing=0.07;
subplot_opts.xmargin=0.04;
subplot_opts.yspacing=0.04;
subplot_opts.ymargin=0.1;
subplot_opts.yoffset=0.03;

tmp1.coeffs=polyfit(tmp1.Ns,mean(tmp1.times*1000,1),1);
disp(tmp1.coeffs);
tmp1.times_fit=polyval(tmp1.coeffs,tmp1.Ns);
subplot_jfm(rr,cc,1,subplot_opts);
plot(tmp1.Ns,mean(tmp1.times*1000,1),'b.',tmp1.Ns,mean(tmp1.times*1000,1),'b',tmp1.Ns,tmp1.times_fit,'r');
xlabel('N = number of samples');
ylabel(ystr);
text(0.05,0.92,'A','FontSize',16,'Units','normalize');
text(0.15,0.7,sprintf('y = %.2f + %.2fx',tmp1.coeffs(2),tmp1.coeffs(1)),'Units','normalize');

tmp2.coeffs=polyfit(tmp2.vals,mean(tmp2.times*1000,1),1);
disp(tmp2.coeffs);
tmp2.times_fit=polyval(tmp2.coeffs,tmp2.vals);
subplot_jfm(rr,cc,2,subplot_opts);
plot(tmp2.vals,mean(tmp2.times*1000,1),'b.',tmp2.vals,mean(tmp2.times*1000,1),'b',tmp2.vals,tmp2.times_fit,'r');
xlabel('K_{true} = true number of clusters');
ylabel(ystr);
text(0.05,0.92,'B','FontSize',16,'Units','normalize');
text(0.15,0.7,sprintf('y = %.2f + %.2fx',tmp2.coeffs(2),tmp2.coeffs(1)),'Units','normalize');

tmp3.coeffs=polyfit(tmp3.vals,mean(tmp3.times*1000,1),1);
disp(tmp3.coeffs);
tmp3.times_fit=polyval(tmp3.coeffs,tmp3.vals);
subplot_jfm(rr,cc,3,subplot_opts);
plot(tmp3.vals,mean(tmp3.times*1000,1),'b.',tmp3.vals,mean(tmp3.times*1000,1),'b',tmp3.vals,tmp3.times_fit,'r');
%plot(tmp3.vals,mean(tmp3.times*1000,1),'b.',tmp3.vals,mean(tmp3.times*1000,1),'b');
xlabel('n = number of dimensions');
ylabel(ystr);
text(0.05,0.92,'C','FontSize',16,'Units','normalize');
text(0.15,0.7,sprintf('y = %.2f + %.2fx',tmp3.coeffs(2),tmp3.coeffs(1)),'Units','normalize');

set(gcf,'paperposition',[0,0,12,3]);
print('../computation_times_01.eps','-depsc2');

end

function simulation_alpha_dependence(simopts)

num_repeats=simopts.num_repeats;

vals=[0.9:0.1:2.3];
mm=length(vals)

if (~exist('simulation_alpha_dependence.mat'))

accs=zeros(1,mm);
times=zeros(1,mm);

for jj=1:mm
	alpha=vals(jj);
	opts=struct;
	opts.ndims=2;
	opts.zdist=1.9;
	opts.num_clusters=6;
	%opts.pops=randi([100,1000],1,opts.num_clusters);
	opts.pops=randi([500,500],1,opts.num_clusters);
	opts.sigma_scale=1;
	opts.spread_factor=0;
	opts.anisotropy_factor=0;
	
	for rr=1:num_repeats
		[clusters,samples,labels]=generate_random_clusters(opts);
		tA=tic;
		labels1=isosplit(samples,struct('K',24,'isocut_threshold',alpha));
		elapsed=toc(tA);
		acc=mean(compute_accuracies(labels,labels1));
		fprintf('alpha=%.3f, %.4f, %.4f s\n',alpha,acc,elapsed);
		times(rr,jj)=elapsed;
		accs(rr,jj)=acc;
	end;
	save('simulation_alpha_dependence.mat','vals','accs','times');
end;
else
	tmp=load('simulation_alpha_dependence.mat');
	vals=tmp.vals;
	accs=tmp.accs;
	times=tmp.times;
end;

disp('==');
accs0=mean(accs,1);
times0=mean(times,1);

for j=1:length(vals)
	fprintf('%.2f & %.0f\\%% & %.2f \\\\\n',vals(j),accs0(j)*100,times0(j));
end;

disp('done with simulation alpha dependence');

end

function simulation_K_dependence(simopts)

num_repeats=simopts.num_repeats;

vals=[3,6,12,24,48,96];
mm=length(vals)

if (~exist('simulation_K_dependence.mat'))

accs=zeros(1,mm);
times=zeros(1,mm);

for jj=1:mm
	opts=struct;
	opts.ndims=2;
	opts.zdist=1.9;
	opts.num_clusters=6;
	%opts.pops=randi([100,1000],1,opts.num_clusters);
	opts.pops=randi([500,500],1,opts.num_clusters);
	opts.sigma_scale=1;
	opts.spread_factor=0;
	opts.anisotropy_factor=0;
	
	for rr=1:num_repeats
		[clusters,samples,labels]=generate_random_clusters(opts);
		tA=tic;
		labels1=isosplit(samples,struct('K',vals(jj)));
		elapsed=toc(tA);
		acc=mean(compute_accuracies(labels,labels1));
		fprintf('K=%.d, %.4f, %.4f s\n',vals(jj),acc,elapsed);
		times(rr,jj)=elapsed;
		accs(rr,jj)=acc;
	end;
	save('simulation_K_dependence.mat','vals','accs','times');
end;
else
	tmp=load('simulation_K_dependence.mat');
	vals=tmp.vals;
	accs=tmp.accs;
	times=tmp.times;
end;

disp('==');
vals;
accs0=mean(accs,1);
times0=mean(times,1);

for j=1:length(vals)
	fprintf('%d & %.0f\\%% & %.2f \\\\\n',vals(j),accs0(j)*100,times0(j));
end;

disp('done with simulation K dependence');

end

function simulation_K_dependence_old

simopts=struct;
simopts.ndims=2;
simopts.zdist=2.5;
simopts.pop_range=[100,1000];
simopts.spread_factor=2;
simopts.anisotropy_factor=1.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simopts.iso_threshold=0.9;
simopts.kmeans_k=[];
simopts.gmm2_k=[];
simopts.hierarchical_cutoff=[];
simopts.db_eps=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simopts.num_trials=100;
simopts.num_tuning_repeats=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simopts.verbose=0;
simopts.seed=1;

simopts.num_clusters=6;

Ks=[3,6,12,24,48,96];

for j=1:length(Ks)
	simopts.K=Ks(j);
	results0=iso_split_simulations(simopts);
	accs(j)=mean(results0{1}.acc);
	elapsed(j)=mean(results0{1}.elapsed);
	disp(accs);
	disp(elapsed);
end;

save('simulation_K_dependence.mat','Ks','accs','elapsed');

disp('==');
Ks
accs
elapsed

disp('done with simulation K dependence');

end

function simulation_m_max_dependence(simopts)

num_repeats=simopts.num_repeats;

vals=[1,5,10,20,50,100,200,500,1000];
mm=length(vals);

if (~exist('simulation_m_max_dependence.mat'))

accs=zeros(1,mm);
times=zeros(1,mm);

for jj=1:mm
	opts=struct;
	opts.ndims=2;
	opts.zdist=1.9;
	opts.num_clusters=6;
	%opts.pops=randi([100,1000],1,opts.num_clusters);
	opts.pops=randi([500,500],1,opts.num_clusters);
	opts.sigma_scale=1;
	opts.spread_factor=0;
	opts.anisotropy_factor=0;
	
	for rr=1:num_repeats
		[clusters,samples,labels]=generate_random_clusters(opts);
		tA=tic;
		labels1=isosplit(samples,struct('K',24,'m_max',vals(jj)));
		elapsed=toc(tA);
		acc=mean(compute_accuracies(labels,labels1));
		fprintf('K=%.d, %.4f, %.4f s\n',vals(jj),acc,elapsed);
		times(rr,jj)=elapsed;
		accs(rr,jj)=acc;
	end;
	save('simulation_m_max_dependence.mat','vals','accs','times');
end;
else
	tmp=load('simulation_m_max_dependence.mat');
	vals=tmp.vals;
	accs=tmp.accs;
	times=tmp.times;
end;

disp('==');
vals;
accs0=mean(accs,1);
times0=mean(times,1);

for j=1:length(vals)
	fprintf('%d & %.0f\\%% & %.2f \\\\\n',vals(j),accs0(j)*100,times0(j));
end;

disp('done with simulation m_max dependence');

end

function simulation_m_max_dependence_old

simopts=struct;
simopts.ndims=2;
simopts.zdist=2.5;
simopts.pop_range=[100,1000];
simopts.spread_factor=2;
simopts.anisotropy_factor=1.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simopts.iso_threshold=0.9;
simopts.kmeans_k=[];
simopts.gmm2_k=[];
simopts.hierarchical_cutoff=[];
simopts.db_eps=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simopts.num_trials=100;
simopts.num_tuning_repeats=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simopts.verbose=0;
simopts.seed=1;

simopts.num_clusters=6;

m_maxs=[1,5,10,20,50,100,200,500,1000];

for j=1:length(m_maxs)
	simopts.m_max=m_maxs(j);
	results0=iso_split_simulations(simopts);
	accs(j)=mean(results0{1}.acc);
	elapsed(j)=mean(results0{1}.elapsed);
	disp(accs);
	disp(elapsed);
end;

save('simulation_m_max_dependence.mat','m_maxs','accs','elapsed');

disp('==');
m_maxs
accs
elapsed

disp('done with simulation m max dependence');

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
