function [algs,algopts]=get_algorithms(gfopts,true_K,dbeps1,dbeps2)

algs={};
algopts={};

if gfopts.do_isosplit
	algs{end+1}=@cluster_alg_iso;
	%algopts{end+1}=struct('name','ISO-SPLIT','K',max(20,true_K*4),'split_threshold',0.9);
    algopts{end+1}=struct('name','ISO-SPLIT','K',max(20,true_K*4),'isocut_threshold',1.2);
end;
if gfopts.do_kmeans
	algs{end+1}=@cluster_alg_kmeans;
	algopts{end+1}=struct('name','K-means (K known)','K',true_K,'r',10);
end;
if gfopts.do_gmm
	algs{end+1}=@cluster_alg_gmm;
	algopts{end+1}=struct('name','Gaussian Mixture (K known)','K',true_K,'r',100,'do_bic_test',0);
end;
if gfopts.do_gmm_bic
	algs{end+1}=@cluster_alg_gmm;
	algopts{end+1}=struct('name','Gaussian Mixture (BIC)','K',true_K,'r',100,'do_bic_test',1);
end;
if (gfopts.do_dbscan)&&(length(dbeps1)>0)
	algs{end+1}=@cluster_alg_dbscan;
	algopts{end+1}=struct('name',sprintf('DBSCAN ($\\epsilon = %g$)',dbeps1),'eps',dbeps1,'minpts',4)
end;
if (gfopts.do_dbscan)&&(length(dbeps2)>0)
	algs{end+1}=@cluster_alg_dbscan;
	algopts{end+1}=struct('name',sprintf('DBSCAN ($\\epsilon = %g$)',dbeps2),'eps',dbeps2,'minpts',4)
end;

end

function labels=cluster_alg_kmeans(samples,opts)
opts.cmethod='k++';
labels=ss_cluster(samples,opts);
end

function labels=cluster_alg_dbscan(samples,opts)
%opts.cmethod='dbd';
opts.cmethod='dbp';
labels=ss_cluster(samples,opts);
end

function labels=cluster_alg_iso(samples,opts)
[labels,info]=isosplit(samples,opts);
%disp(info);
end

function [labels,LL,BIC]=cluster_alg_gmm(samples,opts)

if ((isfield(opts,'do_bic_test'))&&(opts.do_bic_test))
	opts.do_bic_test=0;
	BICs=[];
	kmax=opts.K*2;
	for k=1:kmax
		fprintf('.');
		opts.K=k;
		[~,~,BIC]=cluster_alg_gmm(samples,opts);
		BICs(k)=BIC;
		if ((k>kmax/2)&&(BICs(k)>BICs(k-1))) break; end;
	end
	%fprintf('\n');
	%figure; plot(1:kmax,BICs); drawnow;
	[~,best_k]=min(BICs);
	fprintf(' best_k=%d\n',best_k);
	opts.K=best_k;
	labels=cluster_alg_gmm(samples,opts);
else
	[~,~,priors,LL,posteriors]=vl_gmm(samples,opts.K,'NumRepetitions',20);
	[~,labels]=max(posteriors,[],1);
	dof=10*opts.K-1;
	BIC=-2*LL+dof*(log(size(samples,2)));
end;

end

function [labels,LL,BIC,AIC]=cluster_alg_gmm_matlab(samples,opts)

options=statset;
gm=fitgmdist(samples',opts.K,'Options',options);
labels=cluster(gm,samples');

LL=-gm.NlogL;
BIC=gm.BIC;
AIC=gm.AIC;

end
