function test_2_2_2016

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../franklab/test_hippocampal_01_28_2016/tetrode2_output'];

o_detect.threshold=4;
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=30;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=600;

%% Bandpass filter options
opts_pre.o_filter.samplefreq=30000;
opts_pre.o_filter.freq_min=100;
opts_pre.o_filter.freq_max=10000;
opts_pre.o_filter.outlier_threshold=500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Whitening options
opts_pre.o_whiten=struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre1.mda'],opts_pre.o_filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],opts_pre.o_whiten);
mscmd_detect([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],o_detect);
mscmd_detect([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],[mfile_path,'/tmp_clips.mda'],o_extract_clips);

fprintf('Reading...\n');
%detect=readmda([mfile_path,'/tmp_detect.mda']);
clips=readmda([mfile_path,'/tmp_clips.mda']);
%pre2=readmda([path0,'/pre2.mda']);

[M,T,NC]=size(clips);

fprintf('Computing norms and peaks...\n');
%clip_norms=sqrt(squeeze(sum(sum(clips.^2,1),2)))';
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

all_templates=zeros(M,T,0);
for tt=4:0.5:12
    indsA=find((clip_peaks>tt)&(clip_peaks<tt+0.5));
    fprintf('tt=%g, using %d/%d clips...\n',tt,length(indsA));
    clipsA=clips(:,:,indsA);
    fprintf('Features...\n');
    [FF,subspace]=ms_event_features(clipsA,6);
    fprintf('isosplit...\n');
    labels=isosplit(FF);
    fprintf('K=%d\n',max(labels));
    
    for k=1:max(labels)
        inds_k=find(labels==k);
        FF_k=FF(:,inds_k);
        FF_k_medoid=compute_medoid(FF_k);
        template0=zeros(M,T);
        for jjj=1:size(subspace,3)
            template0=template0+FF_k_medoid(jjj)*subspace(:,:,jjj);
        end;

        all_templates=cat(3,all_templates,template0);
        alpha=0.01;
        inds_outliers=find_outliers(clipsA(:,:,inds_k),template0,alpha);
        fprintf('Found %d/%d outliers\n',length(inds_outliers),length(inds_k));
        labels(inds_k(inds_outliers))=0;
        %FF_tmp=ms_event_features(clipsA(:,:,inds_k),3);
        %figure; ms_view_clusters(FF_tmp,labels(inds_k));
        %figure; ms_view_templates(cat(3,template0,get_example_clips(clipsA(:,:,find(labels==k)),ones(1,length(find(labels==k))),10)));
    end;

end;
all_templates=sort_templates(all_templates);
figure; ms_view_templates(all_templates);

K=size(all_templates,3);
MM=zeros(K,K);
for k1=1:K
    template1=all_templates(:,:,k1);
    for k2=1:K
        template2=all_templates(:,:,k2);
        MM(k1,k2)=compute_correlation(template1(:),template2(:))^2;
    end;
end;

figure; imagesc(MM'); colorbar;
MM2=MM>0.75;
CC=greedy_cliques(MM2);
figure; imagesc(MM2'); colorbar;
% new_templates=zeros(M,T,length(CC));
% for j=1:length(CC)
%     disp(CC{j});
%     new_templates(:,:,j)=mean(all_templates(:,:,CC{j}),3);
% end;
% figure; ms_view_templates(new_templates);

new_templates=zeros(M,T,0);
for j=1:length(CC)
    new_templates=cat(3,new_templates,all_templates(:,:,CC{j}));
    new_templates=cat(3,new_templates,zeros(M,T,2));
end;
figure; ms_view_templates(new_templates);

end

function C=greedy_cliques(A)
for j=1:size(A,1)
    A(j,j)=0;
end;
C={};
while max(A(:))>0
    cc=maximalCliques(A);
    [~,ind0]=max(sum(cc,1));
    ind0=ind0(1);
    C{end+1}=find(cc(:,ind0));
    A(find(cc(:,ind0)),:)=0;
    A(:,find(cc(:,ind0)))=0;
end;
end

function ret=compute_correlation(v1,v2)
ret=sum(v1.*v2)/sqrt(sum(v1.^2)*sum(v2.^2));
end

function outlier_inds=find_outliers(clips,template,alpha)
[M,T,NC]=size(clips);
sumsqr_diffs=reshape(sum(sum((clips-repmat(template,1,1,NC)).^2,1),2),1,NC);
sumsqr_diffs=sumsqr_diffs/(1.2^2);
pp=chi2cdf(sumsqr_diffs,M*T);
outlier_inds=find(pp>1-alpha);
end

function templates=sort_templates(templates)

[M,T,K]=size(templates);

vals=squeeze(sum(templates.^2,2));
norm_vals=sqrt(sum(vals,1));
vals=vals./repmat(norm_vals,M,1);
vals=sum(vals.*repmat((1:M)',1,K),1);

[~,inds]=sort(vals);
templates=templates(:,:,inds);

end

function example_clips=get_example_clips(clips,labels,num_per_label)
[M,T,NC]=size(clips);
K=max(labels);
example_clips=zeros(M,T,K*(num_per_label+1));
for k=1:K
    inds=find(labels==k);
    inds2=inds(randi(length(inds),1,num_per_label));
    i1=(k-1)*(num_per_label+1)+1;
    example_clips(:,:,i1:i1+num_per_label-1)=clips(:,:,inds2);
end;
end


function scores=compute_scores_via_partition(clips)

[M,T,NC]=size(clips);
num_splits=6;
partition=rand_orthant_partition(clips,num_splits);
scores=zeros(1,NC);
L=max(partition);
tA=tic;
for ii=1:L
    if (toc(tA)>1)
        fprintf('ii=%d/%d\n',ii,L);
        tA=tic;
    end;
    inds1=find(partition==ii);
    clips1=clips(:,:,inds1);
    scores(inds1)=compute_scores(clips1);
end;
end

function [scores,clips_norms]=compute_scores(clips,clips_ref)
if nargin<2, clips_ref=clips; end;
[M,T,NC]=size(clips);
[~,~,NC_ref]=size(clips_ref);
ips=zeros(NC_ref,NC);
for j=1:NC_ref
    ips(j,:)=squeeze(sum(sum(clips.*repmat(clips_ref(:,:,j),1,1,NC),1),2));
end;
clips_norms=sqrt(squeeze(sum(sum(clips.^2,1),2)));
clips_ref_norms=sqrt(squeeze(sum(sum(clips_ref.^2,1),2)));
[NORMS1,NORMS2]=ndgrid(clips_ref_norms,clips_norms);
LARGER=max(NORMS1,NORMS2);

scores1=ips./LARGER;
%scores1=ips./(NORMS1.*NORMS2);
%scores1=ips;
scores1_sorted=sort(scores1,1,'descend');
%scores=mean(scores1_sorted(1:20,:),1);
scores=mean(scores1_sorted(20,:),1);

%scores1=(2*ips-NORMS1.^2)/(M*T);
%scores1=ips./(NORMS1.*NORMS2);
%scores1_sorted=sort(scores1,1,'descend');

%scores_a=mean(scores1_sorted(1:20,:),1);
%sigma=sqrt(var(clips(:)));
%scores=2*scores_a/sigma^2;

end

function partition=rand_orthant_partition(clips,num_splits)
[M,T,NC]=size(clips);
partition=zeros(1,NC);

rand_vecs=randn(M,T,num_splits);
ips=zeros(num_splits,NC);
fprintf('Computing ips...\n')
for j=1:num_splits
    ips(j,:)=squeeze(sum(sum(clips.*repmat(rand_vecs(:,:,j),1,1,NC),1),2));
end;

for ii=1:2^num_splits
    to_use=ones(1,NC);
    for j=1:num_splits
        b0=bitget(ii-1,j);
        ip0=ips(j,:);
        if (length(find(to_use==1))>0)
            cutoff=median(ip0(find(to_use==1)));
        else
            cutoff=0;
        end;
        if (b0==0)
            inds0=find(ip0>=cutoff);
        else
            inds0=find(ip0<cutoff);
        end;
        to_use(inds0)=0;
    end;
    partition(find(to_use))=ii;
end;
end

function [clips,inds]=random_orthant(clips,num_splits)

[M,T,NC]=size(clips);

if num_splits==0, return; end;

if num_splits>1
    inds=1:NC;
    for j=1:num_splits
        [clips,inds_tmp]=random_orthant(clips,1);
        inds=inds(inds_tmp);
    end;
    return;
end;

randvec=randn(M,T);
randips=squeeze(sum(sum(repmat(randvec,1,1,NC).*clips,1),2));
[~,sort_inds]=sort(randips,'descend');
inds=sort_inds(1:floor(length(sort_inds)*0.5));
clips=clips(:,:,inds);

end

function m = compute_medoid(X)

[M,N]=size(X);
 
dists=zeros(N,N);
for m=1:M
    [grid1,grid2]=ndgrid(X(m,:),X(m,:));
    dists=dists+(grid1-grid2).^2;
end;
 
avg_dists=mean(dists,1);
[~,ind]=min(avg_dists);
 
m=X(:,ind);
 
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
