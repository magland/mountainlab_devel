function ret=check_gaussian_intersection(mu1,sigma1,mu2,sigma2,rr)

if (nargin<1) test_check_gaussian_intersection; return; end;

% ret=check_gaussian_intersection_old(mu1,sigma1,mu2,sigma2,rr,1000);
% return;

A1=2*inv(sigma1);
b1=-2*inv(sigma1)*mu1;
alpha1=mu1'*inv(sigma1)*mu1-rr^2;

A2=2*inv(sigma2);
b2=-2*inv(sigma2)*mu2;
alpha2=mu2'*inv(sigma2)*mu2-rr^2;

dist0=ellipsoid_distance(A1,b1,alpha1,A2,b2,alpha2);
ret=(dist0==0);

end

function ret=check_gaussian_intersection_old(mu1,sigma1,mu2,sigma2,rr,num_intersection_trials)

if (nargin<1) test_check_gaussian_intersection; return; end;

n=length(mu1);

m=num_intersection_trials;

T1=cholcov(inv(sigma1)); %trans(T1)*T1=inv(sigma1)
T2=cholcov(inv(sigma2)); %trans(T2)*T2=inv(sigma2)

XX=randn(n,m); norms=sqrt(sum(XX.^2,1)); XX=XX./repmat(norms,n,1)*rr;

YY1A=inv(T1)*XX+repmat(mu1,1,m);
YY1B=T2*(YY1A-repmat(mu2,1,m));
YY2A=inv(T2)*XX+repmat(mu2,1,m);
YY2B=T1*(YY2A-repmat(mu1,1,m));

% figure;
% plot_samples(XX,ones(1,size(YY1A,2))*1);
% plot_samples(YY1A,ones(1,size(YY1A,2))*2);
% plot_samples(YY1B,ones(1,size(YY1B,2))*3);
% plot_samples(YY2A,ones(1,size(YY1B,2))*4);
% plot_samples(YY2B,ones(1,size(YY1B,2))*5);
% xlim([-8,8]); ylim([-8,8]);

if (min(sqrt(sum(YY1B.^2,1)),[],2)<=rr)
	ret=1;
	return;
end;
if (min(sqrt(sum(YY2B.^2,1)),[],2)<=rr)
	ret=1;
	return;
end;

ret=0;

end

function samples=generate_samples(mu,sigma,num)

n=length(mu);
samples=randn(n,num);
%norms=sqrt(sum(samples.^2,1));
%samples=samples./repmat(norms,n,1)*1.93;

T=cholcov(inv(sigma));

samples=repmat(mu,1,size(samples,2))+inv(T)*samples;

end

function plot_samples(samples,labels)

colors={'r.','g.','b.','k.','y.','m.','c.','r+','g+','b+','k+','y+','m+','c+'};
for j=1:max(labels)
	xx=samples(1,find(labels==j));
	yy=samples(2,find(labels==j));
	col=colors{mod(j-1,length(colors))+1};
	plot(xx,yy,col);
	if (j==1) hold on; end;
end;
% set(gca,'xtick',[],'ytick',[])
% sub_pos = get(gca,'position'); % get subplot axis position
% sub_pos(1)=sub_pos(1)-sub_pos(3)*0.2;
% sub_pos(3)=sub_pos(3)*1.4;
% sub_pos(2)=sub_pos(2)-sub_pos(4)*0.2;
% sub_pos(4)=sub_pos(4)*1.2;
% set(gca,'position',sub_pos) % stretch its width and height
xmin=min(samples(1,:)); xmax=max(samples(1,:));
ymin=min(samples(2,:)); ymax=max(samples(2,:));
xlim([floor(xmin-2),ceil(xmax+2)]);
ylim([floor(ymin-2),ceil(ymax+2)]);

end

function test_check_gaussian_intersection

close all;

mu1=[0;0]; sigma1=[0.3,0.5;0.5,1];
mu2=[-2;4]; sigma2=[1,-0.5;-0.5,1];

samples1=generate_samples(mu1,sigma1,3000);
samples2=generate_samples(mu2,sigma2,3000);

samples=cat(2,samples1,samples2);
labels=cat(2,ones(1,size(samples1,2))*1,ones(1,size(samples2,2))*2);

plot_samples(samples,labels);

rr=3;
disp(check_gaussian_intersection(mu1,sigma1,mu2,sigma2,rr));

generate_random_clusters;

end
