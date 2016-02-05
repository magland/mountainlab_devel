function tails_toward_zero

%test1;
test2;

end

function test2

close all;

% Set up the true shapes
T=200;
Tcenter=101;
M=2;
tt=linspace(-6,6,T);

pop0=400; pop1=400; pop2=400; pop3=400;
ampl1=8; ampl2=8; ampl3=8;

neurons={};

shape0=zeros(M,T);
shape0(1,:)=ampl1*(abs(tt)<=1).*((1-abs(tt)).^3-0.1);
shape0(2,:)=-ampl1*(abs(tt)<=1)*0.2;
NNN.shape=shape0;
NNN.population=pop1;
neurons{end+1}=NNN;

shape0=zeros(M,T);
shape0(1,:)=ampl2*(abs(tt)<=1)*0.2;
shape0(2,:)=ampl2*(abs(tt)<=1).*((1-abs(tt)).^1);
NNN.shape=shape0;
NNN.population=pop2;
neurons{end+1}=NNN;

shape0=zeros(M,T);
shape0(1,:)=ampl3*(abs(tt)<=1)*0.2;
shape0(2,:)=ampl3*(abs(tt)<=1).*((1-abs(tt)).^3-0.1);
NNN.shape=shape0;
NNN.population=pop3;
neurons{end+1}=NNN;

% Set up the simulated clips and the true labels
clips=zeros(M,T,0);
labels_true=[];

X0=randn(2,T,pop0);
clips=cat(3,clips,X0);
labels_true=[labels_true,ones(1,pop0)*(length(neurons)+1)];

for j=1:length(neurons)
    NNN=neurons{j};
    pop0=NNN.population;
    shape0=NNN.shape;
    shapes=repmat(shape0,1,1,pop0).*repmat(rand(1,1,pop0),M,T,1);
    X1=randn(M,T,pop0)+shapes;
    clips=cat(3,clips,X1);
    labels_true=[labels_true,ones(1,pop0)*j];
end;

N=length(labels_true);

scramble=randsample(N,N);
clips=clips(:,:,scramble);
labels_true=labels_true(scramble);

for j=1:N
    clips(:,:,j)=smooth(clips(:,:,j));
end;

figure;
ms_view_templates(clips(:,:,1:20)); title('Example clips');

figure;
ms_view_templates_from_clips(clips,labels_true); set(gcf,'Name','True templates');

num_features=2;
FF=ms_event_features(clips,num_features);

figure; 
ms_view_clusters(FF,labels_true);
hold on;
th=linspace(0,2*pi,40);
for r0=4:4:20
    plot(r0*cos(th),r0*sin(th),'m--');
    hold on;
end;
title('Clustering in spherical shells');


end

function test1

close all;

% Set up the true shapes
T=200;
Tcenter=101;
M=2;
tt=linspace(-6,6,T);

pop0=400; pop1=400; pop2=400; pop3=400;
ampl1=6; ampl2=6; ampl3=8;

neurons={};

shape0=zeros(M,T);
shape0(1,:)=ampl1*(abs(tt)<=1).*((1-abs(tt)).^3-0.1);
shape0(2,:)=-ampl1*(abs(tt)<=1)*0.2;
NNN.shape=shape0;
NNN.population=pop1;
neurons{end+1}=NNN;

shape0=zeros(M,T);
shape0(1,:)=ampl2*(abs(tt)<=1)*0.2;
shape0(2,:)=ampl2*(abs(tt)<=1).*((1-abs(tt)).^1);
NNN.shape=shape0;
NNN.population=pop2;
neurons{end+1}=NNN;

shape0=zeros(M,T);
shape0(1,:)=ampl3*(abs(tt)<=1)*0.2;
shape0(2,:)=ampl3*(abs(tt)<=1).*((1-abs(tt)).^3-0.1);
NNN.shape=shape0;
NNN.population=pop3;
neurons{end+1}=NNN;

% Set up the simulated clips and the true labels
clips=zeros(M,T,0);
labels_true=[];

X0=randn(2,T,pop0);
clips=cat(3,clips,X0);
labels_true=[labels_true,ones(1,pop0)*(length(neurons)+1)];

for j=1:length(neurons)
    NNN=neurons{j};
    pop0=NNN.population;
    shape0=NNN.shape;
    shapes=repmat(shape0,1,1,pop0).*repmat(rand(1,1,pop0),M,T,1);
    X1=randn(M,T,pop0)+shapes;
    clips=cat(3,clips,X1);
    labels_true=[labels_true,ones(1,pop0)*j];
end;

N=length(labels_true);

scramble=randsample(N,N);
clips=clips(:,:,scramble);
labels_true=labels_true(scramble);

for j=1:N
    clips(:,:,j)=smooth(clips(:,:,j));
end;

figure;
ms_view_templates(clips(:,:,1:20)); title('Example clips');

figure;
ms_view_templates_from_clips(clips,labels_true); set(gcf,'Name','True templates');

num_features=2;
FF=ms_event_features(clips,num_features);

figure; 
ms_view_clusters(FF,labels_true); title('True labels');

labels=isosplit(FF);

figure; 
ms_view_clusters(FF,labels); title('isosplit clustering');

isobranch_opts.isocut_threshold=1.2;
isobranch_opts.K_init=30;
isobranch_opts.min_cluster_split_size=50;
isobranch_opts.num_features=6;
dists=squeeze(max(abs(clips(:,Tcenter,:)),[],1));
labels2=isobranch(clips,dists,isobranch_opts);

figure; 
ms_view_clusters(FF,labels2); title('isobranch clustering');

figure; 
ms_view_templates_from_clips(clips,labels2); set(gcf,'Name','Derived templates');

end

function [labels,clusters]=isobranch(clips,dists,opts)

[M,T,NC]=size(clips);
%Tcenter=ceil((T+1)/2);

clusters={};

isosplit_opts.isocut_threshold=opts.isocut_threshold;
isosplit_opts.K=opts.K_init;

fprintf('features...\n');
features=ms_event_features(clips,opts.num_features);
fprintf('isosplit...\n');
labels=isosplit(features,isosplit_opts);

K=max(labels);

fprintf('%d clusters found (size %d).\n',K,length(labels));

aa=opts.min_cluster_split_size;

for k=1:K
    indices_k=find(labels==k);
    clips0=clips(:,:,indices_k);
    dists0=dists(indices_k);
    if (size(clips0,3)>aa*2)
        [~,sort_inds]=sort(dists0);
        inds_to_use=sort_inds(aa:end);
        clips1=clips0(:,:,inds_to_use);
        [~,clusters2]=isobranch(clips1,dists0(inds_to_use),opts);
        if (length(clusters2)>1)
            for j=1:length(clusters2)
                CC=clusters2{j};
                CC.inds=indices_k(inds_to_use(CC.inds));
                clusters{end+1}=CC;
            end;
        else
            CC.template=mean(clips1,3);
            CC.inds=indices_k;
            clusters{end+1}=CC;
        end;
    else
        CC.template=mean(clips0,3);
        CC.inds=indices_k;
        clusters{end+1}=CC;
    end;
end;

labels=zeros(1,NC);
for j=1:length(clusters)
    CC=clusters{j};
    labels(CC.inds)=j;
end;

end

function inds=biased_sampling(FF,N)
[M,N0]=size(FF);
normsqr=sum(FF.^2,1);

scores=rand(1,N0).*normsqr;
scores_sorted=sort(scores,'descend');
cutoff=scores_sorted(N);

inds=find(scores>=cutoff);
inds=inds(1:N);

end


function inds=biased_sampling_old(FF,N)
[M,N0]=size(FF);
normsqr=sum(FF.^2,1);
[~,sort_inds]=sort(normsqr);
normsqr=normsqr(sort_inds);
CN=cumsum(normsqr);
CN=CN/CN(end);
inds=zeros(1,N);
ii=1;
for j=1:N
    while ((CN(ii)<j/N)&&(ii<N0))
        ii=ii+1;
    end;
    inds(j)=sort_inds(ii);
end;
end

function evt=smooth(evt)
sigma=15;
M=size(evt,1);
T=size(evt,2);
evt_hat=fftshift(fft(fftshift(evt),[],2));
kk=(1:T)-floor((T+1)/2)-1;
evt_hat=evt_hat.*repmat(exp(-kk.^2/(2*sigma^2)),M,1);
evt=fftshift(ifft(fftshift(evt_hat),[],2));
end
