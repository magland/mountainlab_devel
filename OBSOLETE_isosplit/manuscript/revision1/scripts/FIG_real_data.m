function FIG_real_data

addpath('../../../core');

close all;

rng(2);

gfopts.do_isosplit=1;
gfopts.do_kmeans=2;
gfopts.do_gmm=1;
gfopts.do_gmm_bic=1;
gfopts.do_dbscan=1;

disp('loading...');
tmp=load('/home/magland/matlab/scda_ss/data_valid/clips_e_short_th120_3ms_fac3.mat');
X=tmp.X;
X=X(:,:,1:1:end);
disp('extracting features...');
FF=ss_eventfeatures(X);
FF=FF(1:6,:);

fontsize=14;
dbeps1=70;
dbeps2=60;
K0=8;

% FF=FF(:,find(FF(1,:)>-600));
% FF=FF(:,find(FF(3,:)>-600));

viewopts.marker_size=5;
viewopts.create_figure=0;
viewopts.show_legend=0;

camera_position=[-8552,-13978,-8058];
camera_target=[-798,-272,-378];
camera_view_angle=10.15;

N=size(FF,2);
N
ss_view_clusters(FF,ones(1,N),viewopts); drawnow;
set(gca,'CameraPosition',camera_position);
set(gca,'CameraTarget',camera_target);
set(gca,'CameraViewAngle',camera_view_angle);


[algs,algopts]=get_algorithms(gfopts,K0,dbeps1,dbeps2)

minpop=50;

F1=figure('color','w','position',[50,50,1200,800]);
F2=figure('color','w','position',[50,50,1200,800]);

subplot_opts.xmargin=0.04;
subplot_opts.xspacing=0.04;
subplot_opts.ymargin=0.04;
subplot_opts.yspacing=0.06;

for j=1:length(algs)
	alg=algs{j};
	opts=algopts{j}
	str=opts.name;
	opts.name=strrep(opts.name,'K known',sprintf('K = %d',K0));
	str=opts.name;
	str=strrep(str,'$','');
	disp(str);
	tA=tic;
	labels=alg(FF,opts);
	toc(tA)
	labels=order_by_population(labels,minpop);
	if (j==1) labels0=labels;
	else [~,labels]=ss_getbestshuffling(labels,labels0);
	end;
	disp_pops(labels);
	figure(F1); subplot_jfm(2,3,j,subplot_opts);
	ss_view_clusters(FF,labels,viewopts); title(str); drawnow;
	if j==1
		XL=xlim; YL=ylim; ZL=zlim;
	else
		xlim(XL); ylim(YL); zlim(ZL);
	end;
	set(gca,'xtick',[],'ytick',[],'ztick',[]); xlabel(''); ylabel(''); zlabel('');
	set(gca,'CameraPosition',camera_position);
	set(gca,'CameraTarget',camera_target);
	set(gca,'CameraViewAngle',camera_view_angle);
	waveforms=extract_waveforms(X,labels);
	figure(F2); subplot_jfm(2,3,j,subplot_opts);
	ss_view_waveforms(waveforms,viewopts); title(str); drawnow;
end;

set(F1,'paperposition',[0,0,9,5]);
print(F1,'../real_data_1.eps','-depsc2');
set(F2,'paperposition',[0,0,9,5]);
print(F2,'../real_data_2.eps','-depsc2');

end

function labels2=order_by_population(labels,minpop)

K=max(labels);
counts=zeros(1,K);
for j=1:K
	counts(j)=length(find(labels==j));
end;
[~,inds]=sort(counts,'descend');
labels2=zeros(size(labels));
for j=1:K
	if (counts(inds(j))>=minpop)
		labels2(find(labels==inds(j)))=j;
	else
		labels2(find(labels==inds(j)))=0;
	end;
end;

end

function disp_pops(labels)
K=max(labels);
counts=zeros(1,K);
for j=1:K
	ct=length(find(labels==j));
	fprintf('%d:%d  ',j,ct);
end;
fprintf('\n');
end

function W=extract_waveforms(X,labels)
K=max(labels);
[M,T,N]=size(X);
W=zeros(M,T,K);
for j=1:K
	W(:,:,j)=mean(X(:,:,find(labels==j)),3);
end;
end
