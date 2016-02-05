function test_hippocampal_02_04_2016

close all; drawnow;

mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode'];

tetrode_num=1;
path0=[mfile_path,sprintf('/output_tetrode%d',tetrode_num)];
if ~exist(path0,'dir') mkdir(path0); end;
extract_raw_data(raw_path,path0,tetrode_num);

plausibility_threshold=0.6;
merge_threshold=0.8;
tt_list=4:1:12;
dtt=1;
num_features=6;
cross_correlograms_max_dt=6000;
sigma=1.5;

o_mild_filter.samplefreq=30000;
o_mild_filter.freq_min=50;
o_mild_filter.freq_max=10000;
o_mild_filter.outlier_threshold=500;
o_filter.samplefreq=30000;
o_filter.freq_min=100;
o_filter.freq_max=10000;
o_filter.outlier_threshold=500;
o_detect.threshold=tt_list(1);
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=15;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=120;
o_whiten=struct;

mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre1.mda'],o_filter);
mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre0_mild.mda'],o_mild_filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],o_whiten);
mscmd_detect([path0,'/pre2.mda'],[path0,'/detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[path0,'/detect.mda'],[path0,'/clips.mda'],o_extract_clips);

%%KlustaKwik comparison
LL=load([raw_path,'/20151208_NNF_r1_tet16_spikes.mat']);
times=[]; labels=[];
times=[times,LL.clustert_1']; labels=[labels,ones(1,length(LL.clustert_1))*1];
times=[times,LL.clustert_2']; labels=[labels,ones(1,length(LL.clustert_2))*2];
times=[times,LL.clustert_3']; labels=[labels,ones(1,length(LL.clustert_3))*3];
times=[times,LL.clustert_4']; labels=[labels,ones(1,length(LL.clustert_4))*4];
times=[times,LL.clustert_5']; labels=[labels,ones(1,length(LL.clustert_5))*5];
times=[times,LL.clustert_6']; labels=[labels,ones(1,length(LL.clustert_6))*6];
times=[times,LL.clustert_7']; labels=[labels,ones(1,length(LL.clustert_7))*7];
times=(times-LL.starttime)*30000-1e6 +3;
inds000=find((times>0)&(times<25e6)); times=round(times(inds000)); labels=labels(inds000);
[~,sort_inds]=sort(times); times=times(sort_inds); labels=labels(sort_inds);
clusters_kk=zeros(3,length(times)); clusters_kk(2,:)=times; clusters_kk(3,:)=labels;
writemda(clusters_kk,[path0,'/clusters_kk.mda']);
mscmd_templates([path0,'/pre2.mda'],[path0,'/clusters_kk.mda'],[path0,'/templates_kk.mda'],struct('clip_size',200));
mscmd_cross_correlograms([path0,'/clusters_kk.mda'],[path0,'/cross_correlograms_kk.mda'],cross_correlograms_max_dt);
figure; ms_view_templates(readmda([path0,'/templates_kk.mda']));
title('Templates from KlustaKwik'); drawnow;
pre2=readmda([path0,'/pre2.mda']);
clips_kk=ms_extract_clips(pre2,times,o_extract_clips.clip_size);
[clips0_kk,clips0_index_kk]=ms_create_clips_index(clips_kk,labels);
writemda(clips0_kk,[path0,'/clips_kk.mda']);
writemda(clips0_index_kk,[path0,'/clips_index_kk.mda']);
view_params.raw=[path0,'/pre2.mda'];
view_params.clusters=[path0,'/clusters_kk.mda'];
view_params.cross_correlograms=[path0,'/cross_correlograms_kk.mda'];
view_params.templates=[path0,'/templates_kk.mda'];
view_params.clips=[path0,'/clips_kk.mda'];
view_params.clips_index=[path0,'/clips_index_kk.mda'];
ms_mountainview(view_params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre2=readmda([path0,'/pre2.mda']);
%spikespy({pre2,times,labels});

fprintf('Reading...\n');
clips=readmda([path0,'/clips.mda']);
[M,T,NC]=size(clips);

fprintf('Shell cluster...\n');
all_templates=shell_cluster(clips,tt_list,dtt,num_features);
all_templates=group_templates(all_templates);
figure; ms_view_templates(all_templates); drawnow;
[~,~,K_all]=size(all_templates);
title('All templates');

fprintf('Posterior probabilities...');
[probs,plausibility_factors]=compute_posterior_probabilities(clips,all_templates,sigma);
corr_matrix=corrcoef(probs');
figure; imagesc(corr_matrix');
title('Correlation matrix for posterior classification prob vectors'); drawnow;
plausibility_factors=probs./repmat(max(probs,[],2),1,size(probs,2));
figure; hist(max(plausibility_factors,[],1),linspace(-0.02,1.02,1000));
title('Plausibility factors for all events'); drawnow;

fprintf('Computing times/labels and coincidence matrix...\n');
detect=readmda([path0,'/detect.mda']);
times0=detect(2,:);
times=[]; labels=[]; used=zeros(size(times0));
%plausible_events=(plausibility_factors>plausibility_threshold);

for k=1:size(plausibility_factors,1);
    inds_k=find((plausibility_factors(k,:)>=plausibility_threshold)&(probs(k,:)>=1.0*max(probs,[],1)));
    times=[times,times0(inds_k)];
    labels=[labels,ones(1,length(inds_k))*k];
    used(inds_k)=1;
end;
[~,sort_inds]=sort(times); times=times(sort_inds); labels=labels(sort_inds);
CM=compute_coincidence_matrix(times,labels);
figure; imagesc(CM'); title('Coincidence matrix'); drawnow;

% AA=all_templates;
% AA=AA-repmat(mean(AA,2),1,T,1);
% AA=reshape(AA,M*T,K_all);
% cmat=corrcoef(AA);
% figure; imagesc(cmat');
% title('Correlation matrix for all templates');
% CC=greedy_cliques(cmat>=merge_threshold);

%CC=greedy_cliques((CM+CM')/2>=merge_threshold);
%CC=greedy_cliques(max(CM,CM')>=merge_threshold);
% 
% times=[];
% labels=[];
% for j=1:length(CC)
%     inds00=CC{j};
%     inds1=find(max(plausibility_factors(inds00,:),[],1)>=plausibility_threshold);
%     times=[times,times0(inds1)];
%     labels=[labels,ones(1,length(inds1))*j];
% end;
% [~,sort_inds]=sort(times); times=times(sort_inds); labels=labels(sort_inds);
% CM=compute_coincidence_matrix(times,labels);
% figure; imagesc(CM'); title('Coincidence matrix'); drawnow;

clusters=zeros(3,length(times));
clusters(2,:)=times;
clusters(3,:)=labels;

fprintf('Writing output and preparing view...\n');
pre2=readmda([path0,'/pre2.mda']);
[clips1,clips1_index]=ms_create_clips_index(ms_extract_clips(pre2,times,o_extract_clips.clip_size),labels);
writemda(clips1,[path0,'/clips.mda']);
writemda(clips1_index,[path0,'/clips_index.mda']);
writemda(clusters,[path0,'/clusters.mda']);

mscmd_cross_correlograms([path0,'/clusters.mda'],[path0,'/cross_correlograms.mda'],cross_correlograms_max_dt);
mscmd_templates([path0,'/pre0_mild.mda'],[path0,'/clusters.mda'],[path0,'/templates_raw.mda'],struct('clip_size',200));
mscmd_templates([path0,'/pre2.mda'],[path0,'/clusters.mda'],[path0,'/templates.mda'],struct('clip_size',200));

templates_raw=readmda([path0,'/templates_raw.mda']);
figure; ms_view_templates(templates_raw);

view_params.raw=[path0,'/pre2.mda'];
view_params.clusters=[path0,'/clusters.mda'];
view_params.cross_correlograms=[path0,'/cross_correlograms.mda'];
view_params.templates=[path0,'/templates.mda'];
view_params.clips=[path0,'/clips.mda'];
view_params.clips_index=[path0,'/clips_index.mda'];
ms_mountainview(view_params);

mscmd_confusion_matrix([path0,'/clusters_kk.mda'],[path0,'/clusters.mda'],[path0,'/confusion_matrix.mda'],4);
CM=readmda([path0,'/confusion_matrix.mda']);
CM1=CM./repmat(sum(CM,1),size(CM,1),1);
CM2=CM./repmat(sum(CM,2),1,size(CM,2));
figure; imagesc(CM1'); colorbar; xlabel('KlustaKwik Clusters'); ylabel('MountainSort Clusters'); title('Confusion Matrix Normalized by Rows');
figure; imagesc(CM2'); colorbar; xlabel('KlustaKwik Clusters'); ylabel('MountainSort Clusters'); title('Confusion Matrix Normalized by Columns');

end

function templates=merge_templates(templates,AM)
[M,T,K]=size(templates);
cliques=greedy_cliques(AM);
templates_new=zeros(M,T,length(cliques));
for cc=1:length(cliques)
    inds00=cliques{cc};
    templates_new(:,:,cc)=mean(templates(:,:,inds00),3);
end;
templates=templates_new;
end

function [CM_normalized,CM]=compute_coincidence_matrix(times,labels)
K=max(labels);
N=max(times);
CM=zeros(K,K);
for k1=1:K
    for k2=1:K
        CM(k1,k2)=length(intersect(times(find(labels==k1)),times(find(labels==k2))));
    end;
end;
%Now, normalize
for k1=1:K
    for k2=1:K
        CM_normalized(k1,k2)=CM(k1,k2)/CM(k1,k1);
    end;
end;
end

function [posterior_probs,corr_matrix]=compute_posterior_probabilities_try(clips,templates,sigma)
fprintf('Computing ips... ');
[M,T,K]=size(templates);
[~,~,NC]=size(clips);
ips=zeros(K,NC);
for k=1:K
    fprintf('%d ',k);
    ips(k,:)=squeeze(sum(sum(clips.*repmat(templates(:,:,k),1,1,NC),1),2));
end;
fprintf('\nComputing norms, etc...\n');
template_norms=reshape(sqrt(sum(sum(templates.^2,1),2)),1,K);
clip_norms=reshape(sqrt(sum(sum(clips.^2,1),2)),1,NC);
diffsqr=repmat(clip_norms,K,1).^2-2*ips+repmat(template_norms',1,NC).^2;

logprobs1=-2*diffsqr/sigma^2/(M*T); %KxNC
logprobs2=-2*clip_norms.^2/sigma^2/(M*T); %1xNC
logprobs3=log(0.05)*ones(1,NC); %1xNC

posterior_probs=zeros(K,NC);
for k=1:K
    posterior_probs(k,:)=exp(logprobs1(k,:) - logsumexp(cat(1,logprobs1(k,:),logprobs2,logprobs3),1) );
end;

corr_matrix=corrcoef(posterior_probs');
end

function [probs,corr_matrix]=compute_posterior_probabilities_with_scaling(clips,templates,sigma)
fprintf('Computing ips... ');
[M,T,K]=size(templates);
[~,~,NC]=size(clips);
ips=zeros(K,NC);
for k=1:K
    fprintf('%d ',k);
    ips(k,:)=squeeze(sum(sum(clips.*repmat(templates(:,:,k),1,1,NC),1),2));
end;
fprintf('\nComputing norms, etc...\n');
template_norms=reshape(sqrt(sum(sum(templates.^2,1),2)),1,K);
clip_norms=reshape(sqrt(sum(sum(clips.^2,1),2)),1,NC);

%scale_constants=ips./repmat((template_norms').^2,1,NC);
%adjusted_template_norms=scale_constants.*repmat((template_norms'),1,NC);
diffsqr=repmat(clip_norms,K,1).^2-ips.^2./repmat((template_norms'),1,NC).^2;

logprobs1=-2*diffsqr/sigma^2/(M*T);
logmargins=logsumexp(cat(1,logprobs1,log(0.2*ones(1,NC))),1);
%logmargins=logsumexp(logprobs1);
logprobs=logprobs1-repmat(logmargins,K,1);
probs=exp(logprobs);
corr_matrix=corrcoef(probs');
end

function [probs,plausibility_factors]=compute_posterior_probabilities(clips,templates,sigma)
fprintf('Computing ips... ');
[M,T,K]=size(templates);
[~,~,NC]=size(clips);
ips=zeros(K,NC);
for k=1:K
    fprintf('%d ',k);
    ips(k,:)=squeeze(sum(sum(clips.*repmat(templates(:,:,k),1,1,NC),1),2));
end;
fprintf('\nComputing norms, etc...\n');
template_norms=reshape(sqrt(sum(sum(templates.^2,1),2)),1,K);
clip_norms=reshape(sqrt(sum(sum(clips.^2,1),2)),1,NC);
diffsqr=repmat(clip_norms,K,1).^2-2*ips+repmat(template_norms',1,NC).^2;

logprobs1=-2*diffsqr/sigma^2/(M*T);
logmargins=logsumexp(cat(1,logprobs1,log(0.2*ones(1,NC))),1);
%logmargins=logsumexp(logprobs1);
logprobs=logprobs1-repmat(logmargins,K,1);
probs=exp(logprobs);

normalizations=exp(logsumexp(cat(1,zeros(1,K),log(0.2*ones(1,K)))));
plausibility_factors=probs./repmat(normalizations',1,NC);

end

function templates=shell_cluster(clips,tt_list,dtt,num_features)
[M,T,NC]=size(clips);

fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

clusterings={};
for ii=1:length(tt_list)
    tt=tt_list(ii);
    tt2=tt+dtt;
    if (ii==length(tt_list)), tt2=inf; end;
    fprintf('tt=%g... ',tt);
    inds_tt=find((clip_peaks>=tt)&(clip_peaks<=tt2));
    CC.tt=tt;
    CC.inds=inds_tt;
    if (length(inds_tt)>1)
        clips_tt=clips(:,:,inds_tt);
        fprintf('features... ');
        [FF_tt,subspace_tt]=ms_event_features(clips_tt,num_features);
        fprintf('isosplit... ');
        labels_tt=isosplit(FF_tt);
        K=max(labels_tt);
        fprintf('K=%d\n',K);
        CC.labels=labels_tt;
        CC.features=FF_tt;
        CC.clips=clips_tt;
        CC.medoids=zeros(M,T,K);
        for k=1:K
            FF_tt_k=FF_tt(:,find(labels_tt==k));
            FF_tt_k_medoid=compute_medoid(FF_tt_k);
            template0=zeros(M,T);
            for jjj=1:size(subspace_tt,3)
                template0=template0+FF_tt_k_medoid(jjj)*subspace_tt(:,:,jjj);
            end;
            CC.medoids(:,:,k)=template0;
        end;
    else
        CC.labels=[];
        CC.features=zeros(num_features,0);
        CC.clips=zeros(M,T,0);
        CC.medoids=zeros(M,T,0);
    end;
    clusterings{end+1}=CC;
end

templates=zeros(M,T,0);
for ii=1:length(clusterings)
    templates=cat(3,templates,clusterings{ii}.medoids);
end;

end

function templates=group_templates(templates)

[M,T,K]=size(templates);
peaks_pos=squeeze(max(templates(:,T/2+1,:),[],1))';
peaks_neg=-squeeze(max(-templates(:,T/2+1,:),[],1))';
peaks=peaks_pos.*(abs(peaks_pos)>abs(peaks_neg))+peaks_neg.*(abs(peaks_pos)<abs(peaks_neg));

[~,sort_inds]=sort(peaks);
templates=templates(:,:,sort_inds);

[M,T,K]=size(templates);
MM=zeros(K,K);
templates=templates-repmat(mean(templates,2),1,T,1);
for k1=1:K
    template1=templates(:,:,k1);
    for k2=1:K
        template2=templates(:,:,k2);
        MM(k1,k2)=compute_correlation(template1(:),template2(:));
    end;
end;

%figure; imagesc(MM'); colorbar;
MM2=MM>0.8;
CC=greedy_cliques(MM2);
%figure; imagesc(MM2'); colorbar;

new_templates=zeros(M,T,0);
for j=1:length(CC)
    new_templates=cat(3,new_templates,templates(:,:,CC{j}));
    %new_templates=cat(3,new_templates,zeros(M,T,2));
end;
templates=new_templates;

end


function extract_raw_data(raw_path,output_path,tetrode_num)

raw_mat_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mat',raw_path);
raw_mda_fname=sprintf('%s/dl12_20151208_NNF_r1_tet16_17.mda',raw_path);
tetrode_fname=sprintf('%s/pre0.mda',output_path);

if (~exist(raw_mda_fname,'file'))
    fprintf('Loading raw data...\n');
    L=load(raw_mat_fname);
    raw=L.dl12_20151208_NNF_r1_tet16_17.channelData';
    fprintf('Writing raw data...\n');
    writemda(raw,raw_mda_fname);
end;

if (~exist(tetrode_fname,'file'))
    fprintf('Reading raw data...\n');
    raw=readmda(raw_mda_fname);

    if (tetrode_num==1)
        tetrode=raw([1,3:6],(1e6+1):26e6);
        tetrode=tetrode(2:end,:)-repmat(tetrode(1,:),size(tetrode,1)-1,1);
    elseif (tetrode_num==2)
        tetrode=raw([1,7:10],(1e6+1):26e6);
        tetrode=tetrode(2:end,:)-repmat(tetrode(1,:),size(tetrode,1)-1,1);
    end;
    fprintf('Writing tetrode data...\n');
    writemda(tetrode,tetrode_fname);
end;

L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];
writemda(L,[output_path,'/locations.mda']);


end

function m = compute_medoid(X)
[M,N]=size(X);
dists=zeros(N,N);
for m=1:M
    [grid1,grid2]=ndgrid(X(m,:),X(m,:));
    dists=dists+sqrt((grid1-grid2).^2);
    %dists=dists+(grid1-grid2).^2;
end;
avg_dists=mean(dists,1);
[~,ind]=min(avg_dists);
m=X(:,ind); 
end


function ret=compute_correlation(v1,v2)
ret=sum(v1.*v2)/sqrt(sum(v1.^2)*sum(v2.^2));
end


function C=greedy_cliques(A)
C={};
used=zeros(1,size(A,1));
while max(A(:))>0
    A0=A; for j=1:size(A,1) A0(j,j)=0; end;
    cc=maximalCliques(A0);
    [maxval,ind0]=max(sum(cc,1));
    if (maxval==1)
        inds0=find(used==0);
        for aa=1:length(inds0)
            C{end+1}=[inds0(aa)];
        end;
        break;
    end;
    ind0=ind0(1);
    ind_cc=find(cc(:,ind0));
    C{end+1}=ind_cc;
    A(ind_cc,:)=0;
    A(:,ind_cc)=0;
    used(ind_cc)=1;
end;
minvals=[];
for j=1:length(C)
    minvals=[minvals,min(C{j})];
end;
[~,sort_inds]=sort(minvals);
C=C(sort_inds);
end

function [ MC ] = maximalCliques( A, v_str )
%MAXIMALCLIQUES Find maximal cliques using the Bron-Kerbosch algorithm
%   Given a graph's boolean adjacency matrix, A, find all maximal cliques 
%   on A using the Bron-Kerbosch algorithm in a recursive manner.  The 
%   graph is required to be undirected and must contain no self-edges.
%
%   V_STR is an optional input string with the version of the Bron-Kerbosch 
%   algorithm to be used (either 'v1' or 'v2').  Version 2 is faster (and 
%   default), and version 1 is included for posterity.
%
%   MC is the output matrix that contains the maximal cliques in its 
%   columns.
%
%   Note: This function can be used to compute the maximal independent sets
%   of a graph A by providing the complement of A as the input graph.  
%
%   Note: This function can be used to compute the maximal matchings of a 
%   graph A by providing the complement of the line graph of A as the input
%   graph.
%
%   Ref: Bron, Coen and Kerbosch, Joep, "Algorithm 457: finding all cliques
%   of an undirected graph", Communications of the ACM, vol. 16, no. 9, 
%   pp: 575–577, September 1973.
%
%   Ref: Cazals, F. and Karande, C., "A note on the problem of reporting 
%   maximal cliques", Theoretical Computer Science (Elsevier), vol. 407,
%   no. 1-3, pp: 564-568, November 2008.
%
%   Jeffrey Wildman (c) 2011
%   jeffrey.wildman@gmail.com
%   
%   Updated: 10/27/2011 - updated documentation & removal of ~ punctuation 
%   to ignore function output arguments for better compatibility with older
%   MATLAB versions prior to 2009b (Thanks to Akli Benali).


% first, some input checking

if size(A,1) ~= size(A,2)
    error('MATLAB:maximalCliques', 'Adjacency matrix is not square.');
elseif ~all(all((A==1) | (A==0)))
    error('MATLAB:maximalCliques', 'Adjacency matrix is not boolean (zero-one valued).')
elseif ~all(all(A==A.'))
    error('MATLAB:maximalCliques', 'Adjacency matrix is not undirected (symmetric).')
elseif trace(abs(A)) ~= 0
    error('MATLAB:maximalCliques', 'Adjacency matrix contains self-edges (check your diagonal).');
end
    
if ~exist('v_str','var')
    v_str = 'v2';
end

if ~strcmp(v_str,'v1') && ~strcmp(v_str,'v2')
    warning('MATLAB:maximalCliques', 'Version not recognized, defaulting to v2.');
    v_str = 'v2';
end


% second, set up some variables

n = size(A,2);      % number of vertices
MC = [];            % storage for maximal cliques
R = [];             % currently growing clique
P = 1:n;            % prospective nodes connected to all nodes in R
X = [];             % nodes already processed


% third, run the algorithm!

if strcmp(v_str,'v1')
    BKv1(R,P,X);
else
    BKv2(R,P,X);
end
    

    % version 1 of the Bron-Kerbosch algo 
    function [] = BKv1 ( R, P, X )
        
        if isempty(P) && isempty(X)
            % report R as a maximal clique
            newMC = zeros(1,n);
            newMC(R) = 1;                   % newMC contains ones at indices equal to the values in R   
            MC = [MC newMC.'];
        else
            for u = P
                P = setxor(P,u);
                Rnew = [R u];
                Nu = find(A(u,:));
                Pnew = intersect(P,Nu);
                Xnew = intersect(X,Nu);
                BKv1(Rnew, Pnew, Xnew);
                X = [X u];
            end
        end
        
    end % BKv1


	% version 2 of the Bron-Kerbosch algo
    function [] = BKv2 ( R, P, X )

        ignore = [];                        % less elegant ignore function output variable, but works with older versions of MATLAB: <2009b
        if (isempty(P) && isempty(X))
            % report R as a maximal clique
            newMC = zeros(1,n);
            newMC(R) = 1;                   % newMC contains ones at indices equal to the values in R   
            MC = [MC newMC.'];
        else
            % choose pivot
            ppivots = union(P,X);           % potential pivots
            binP = zeros(1,n);
            binP(P) = 1;                    % binP contains ones at indices equal to the values in P          
            % rows of A(ppivots,:) contain ones at the neighbors of ppivots
            pcounts = A(ppivots,:)*binP.';  % cardinalities of the sets of neighbors of each ppivots intersected with P
            [ignore,ind] = max(pcounts);
            u_p = ppivots(ind);             % select one of the ppivots with the largest count
            
            for u = intersect(find(~A(u_p,:)),P)   % all prospective nodes who are not neighbors of the pivot
                P = setxor(P,u);
                Rnew = [R u];
                Nu = find(A(u,:));
                Pnew = intersect(P,Nu);
                Xnew = intersect(X,Nu);
                BKv2(Rnew, Pnew, Xnew);
                X = [X u];
            end
        end
        
    end % BKv2
       

end % maximalCliques
