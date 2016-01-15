function test_004

close all;

mfile_path=fileparts(mfilename('fullpath'));

ch0=1;
AM=readmda(sprintf('%s/../example_data/adjacency.mda',mfile_path));
ch=find(AM(ch0,:));
disp(ch);

fprintf('Reading X...\n');
X=readmda(sprintf('%s/../example_data/filt2_white.mda',mfile_path));
X=X(ch,:);
[M,N]=size(X);

opts.detect_interval=200;
opts.detect_threshold=5;
opts.clip_size=50;

X0=X;

fprintf('Detecting events...\n');
times=detect_events(X0,find(ch==ch0),opts);
fprintf('Detected %d events.\n',length(times));

fprintf('Extract clips...\n');
clips=ms_extract_clips(X,times,opts.clip_size);
ss_view_waveforms(clips(:,:,1:40));
drawnow;

fprintf('Finding non-zero events...\n');
npca=10;
num_neighbor_pts=25;
iii=find_nonzero_events(clips,npca,num_neighbor_pts);
clips=clips(:,:,iii);
times=times(iii);
fprintf('Using %d non-zero events...\n',length(iii));
ss_view_waveforms(clips(:,:,1:40));
drawnow;

fprintf('ms_event_features...\n');
FF=ms_event_features(clips,6);
fprintf('isosplit...\n');
labels=isosplit(FF);

ss_view_clusters(FF,labels);
drawnow;

K=max(labels);
templates=zeros(M,opts.clip_size,K);
for k=1:K
    clips0=clips(:,:,find(labels==k));
    templates(:,:,k)=mean(clips0,3);
end;

ss_view_waveforms(templates);
drawnow;

spikespy({X,times,labels});

end

function X=zero_out_events(X,times,clip_size)

[M,N]=size(X);
T=clip_size;
C=length(times);

tt1=-ceil((clip_size)/2);
tt2=tt1+clip_size-1;
if (min(times+tt1)<1) error('Invalid time in zero_out_events'); end;
if (max(times+tt2)>N) error('Invalid time in zero_out_events'); end;
for j=1:C
	X(:,times(j)+tt1:times(j)+tt2)=0;
end;

end

function times=detect_events(X,ch,opts)

fprintf('ms_detect...\n');
times=ms_detect(X(ch,:),opts);

end

function inds=find_nonzero_events(clips,npca,num_neighbor_pts)
V=reshape(clips,size(clips,1)*size(clips,2),size(clips,3));
fprintf('PCA...\n');
FF=pca_features(V);
FF=FF(1:npca,:);

fprintf('find_nonzero_vectors...\n');
inds=find_nonzero_vectors(FF,num_neighbor_pts);

end

function inds=find_nonzero_vectors(V,num_neighbor_pts)

cutoff=3;

[M,N]=size(V);
sigma_guess=sqrt(var(V(:)));
norms=sqrt(sum(V.^2,1));

fprintf('knnsearch...\n');
[ids,ds]=knnsearch(V',V','K',num_neighbor_pts);

% Volume of sphere is (pi)^[M/2] / gamma((M+2)/2) * R^M
% where R=ds(:,end)
log_est_densities=log(num_neighbor_pts)-(M/2)*log(pi)+gammaln((M+2)/2)-M*log(ds(:,end)');

fprintf('estimate_sigma...\n');
[est_sigma,est_log_density_at_zero]=estimate_sigma(M,norms,log_est_densities,sigma_guess);

% PDF of gaussian: (2*pi)^(-M/2)*sigma^M*exp(-norm^2/(2*sigma^2))
%log_expected_densities=log(N)-M/2*log(2*pi)-M*log(est_sigma)-1/2*norms.^2/est_sigma^2;
log_expected_densities=est_log_density_at_zero-1/2*(norms).^2/est_sigma^2;

figure; plot(log_expected_densities,log_est_densities,'b.');
xlabel('Log expected densities');
ylabel('Log estimated densities');
hold on;
plot(xlim,xlim,'r');

scores=log_est_densities-log_expected_densities;
figure;
plot(norms,scores,'b.');
hold on;
plot(norms(scores>cutoff),scores(scores>cutoff),'r.');
xlabel('Norm');
ylabel('Score');

figure; set(gcf,'position',[100,1200,1000,500]);
subplot(1,2,1);
bins=linspace(min(scores),max(scores),1000);
hist(scores,bins); hold on;
hist(scores(scores>cutoff),bins); hold on;
xlabel('Scores');

subplot(1,2,2);
bins=linspace(min(norms),max(norms),1000);
hist(norms,bins); hold on;
hist(norms(scores>cutoff),bins); hold on;
xlabel('Norms');

inds=find(scores>cutoff);

end

function [est_sigma,est_log_density_at_zero]=estimate_sigma(M,norms,log_est_densities,sigma_guess)
%figure; plot(norms.^2,log_est_densities,'b.');
%xlabel('norms^2');
%ylabel('log_est_densities');
inds=find(norms<=sigma_guess*sqrt(M));
[slope,intercept]=estimate_slope_intercept(norms(inds).^2,log_est_densities(inds));
est_sigma=sqrt(-1/(2*slope));
est_log_density_at_zero=intercept;
end

function [slope,intercept]=estimate_slope_intercept(X,Y)
coeffs=polyfit(X,Y,1);
slope=coeffs(1);
intercept=coeffs(2);
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