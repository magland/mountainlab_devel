function isosplit_tests

close all;

run_test_A(0);
run_test_A(4);
run_test_A(8);
run_test_A(12);

end

function run_test_A(x_offset)

seed0=randi(10000);
seed0=6245;
rng(seed0);
centers={[0,0],[x_offset,5]};
pops={1000,1000};
shapes={[3,1,0],[3,1,0]};
opts.K=10;
opts.use_geometric_median=0;
opts.verbose2=1;

fprintf('seed = %d\n',seed0);

samples=zeros(2,0);
true_labels=zeros(1,0);

for j=1:length(centers)
	xx=randn(1,pops{j});
	yy=randn(1,pops{j});
	shape=shapes{j};
	xx2=xx*shape(1)+yy*shape(3);
	yy2=yy*shape(2)-xx*shape(3);
	center=centers{j};
	xx2=xx2+center(1);
	yy2=yy2+center(2);
	tmp=zeros(2,pops{j});
	tmp(1,:)=xx2; tmp(2,:)=yy2;
	samples=[samples,tmp];
	true_labels=[true_labels,ones(1,pops{j})*j];
end;

colors='rgbkymc';
[labels,info]=isosplit_dev(samples,opts);
fprintf('num clusters = %d\n',max(labels));


end

function test_split_technique(samples)

[M,N]=size(samples);

AA=zeros(5,N);
AA(1:2,:)=samples;
AA(3,:)=samples(1,:).*samples(1,:);
AA(4,:)=samples(1,:).*samples(2,:);
AA(5,:)=samples(2,:).*samples(2,:);

norms=sqrt(sum(AA.^2,2));
AA=AA./repmat(norms,1,N);

FF=ms_event_features(reshape(AA,size(AA,1),1,size(AA,2)),3);
labels=isosplit_dev(FF);
figure; ms_view_clusters(samples,labels); axis square; daspect([1,1,1]);

end

