function test_003

close all;

mfile_path=fileparts(mfilename('fullpath'));

ch0=12;
AM=readmda(sprintf('%s/../example_data/adjacency.mda',mfile_path));
ch=find(AM(ch0,:));
disp(ch);

X=readmda(sprintf('%s/../example_data/filt2_white.mda',mfile_path));
X=X(ch,:);
[M,N]=size(X);

if 1
    fprintf('Extract clips...\n');
    opts.detect_interval=500;
    opts.detect_threshold=4;
    opts.clip_size=50;
    fprintf('ms_detect...\n');
    [Tpos,Tneg]=ms_detect(X(find(ch==ch0),:),opts);
    times=sort([Tpos,Tneg]);
    fprintf('ms_extract_clips...\n');
    clips=ms_extract_clips(X,times,opts.clip_size);
    writemda(clips,'scratch/clips.mda');
end;
fprintf('Reading clips...\n');
clips=readmda('scratch/clips.mda');
[M,T,NC]=size(clips);
ss_view_waveforms(clips(:,:,1:20));
fprintf('%d clips...\n',size(clips,3));

V=reshape(clips,size(clips,1)*size(clips,2),size(clips,3));
fprintf('PCA...\n');
FF=pca_features(V);
FF=FF(1:6,:);
disp('knnsearch...');
[idx,d]=knnsearch(FF',FF','K',200);
tmp=zeros(size(V));
dists=[];
sum_dists=zeros(1,NC);
num_dists=0;
for aa=2:size(idx,2)
    diff0=V(:,idx(:,aa))-V;
    dists0=sqrt(sum(diff0.^2,1));
    sum_dists=sum_dists+dists0;
    num_dists=num_dists+1;
end;
scores=sqrt(sum(V.^2,1))-sum_dists/num_dists;

cutoff=find_cutoff(scores,4);
%cutoff=4000;
figure; hist(scores,1000); hold on;
plot([cutoff,cutoff],ylim);

iii=find(scores>cutoff);
ss_view_waveforms(clips(:,:,iii(1:40)));

W=V(:,iii);
FF2=pca_features(W);
FF2=FF2(1:3,:);
labels=isosplit(FF2);
ss_view_clusters(FF2,labels);

K=max(labels);
templates=zeros(M,T,K);
for k=1:K
    clips0=clips(:,:,iii(find(labels==k)));
    templates(:,:,k)=mean(clips0,3);
end;

ss_view_waveforms(templates);

end

function cutoff=find_cutoff(scores,num_stdev)
[mu_est, sigma_est, p_est]=gaussian_mixture_model(scores,3);
[~,ind]=max(p_est);
cutoff=mu_est(ind)+sigma_est(ind)*num_stdev;
return;

for pass=1:5
    mu=mean(scores);
    sigma=sqrt(var(scores));
    cutoff=mu+sigma*num_stdev;
    scores=scores(abs(scores)<cutoff);
end;
end

function FF=pca_features(X)
[U,D]=eig(X*X');
[d,I]=sort(diag(D),'descend');
U=U(:,I);
FF=U'*X;
end

