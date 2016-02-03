function test_2_2_2016_B

mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../franklab/test_hippocampal_01_28_2016/tetrode2_output'];

o_filter.samplefreq=30000;
o_filter.freq_min=100;
o_filter.freq_max=10000;
o_filter.outlier_threshold=500;
o_detect.threshold=4;
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=30;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=60;
o_whiten=struct;

mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre1.mda'],o_filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],o_whiten);
mscmd_detect([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],[mfile_path,'/tmp_clips.mda'],o_extract_clips);

fprintf('Reading...\n');
clips=readmda([mfile_path,'/tmp_clips.mda']);

[M,T,NC]=size(clips);

fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

tt_list=4:0.5:12;
dt=0.5;
num_features=6;

clusterings={};
for ii=1:length(tt_list)
    tt=tt_list(ii);
    tt2=tt+dt;
    if (ii==length(tt_list)) tt2=inf; end;
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

all_templates=zeros(M,T,0);
for ii=1:length(clusterings)
    all_templates=cat(3,all_templates,clusterings{ii}.medoids);
end;

all_templates=group_templates(all_templates);

figure; ms_view_templates(all_templates);
K=size(all_templates,3);

fprintf('Computing ips... ');
ips=zeros(K,NC);
for k=1:K
    fprintf('%d ',k);
    ips(k,:)=squeeze(sum(sum(clips.*repmat(all_templates(:,:,k),1,1,NC),1),2));
end;
fprintf('\nComputing norms, etc...\n');
template_norms=reshape(sqrt(sum(sum(all_templates.^2,1),2)),1,K);
clip_norms=reshape(sqrt(sum(sum(clips.^2,1),2)),1,NC);
diffsqr=repmat(clip_norms,K,1).^2-2*ips+repmat(template_norms',1,NC).^2;

sigma=1.2;
logprobs=-2*diffsqr/sigma^2;
logmargins=logsumexp(cat(1,logprobs,-10*ones(1,NC)),1);
logprobs=logprobs-repmat(logmargins,K,1);
probs=exp(logprobs);

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

function templates=group_templates(templates)

[M,T,K]=size(templates);
MM=zeros(K,K);
for k1=1:K
    template1=templates(:,:,k1);
    for k2=1:K
        template2=templates(:,:,k2);
        MM(k1,k2)=compute_correlation(template1(:),template2(:));
    end;
end;

figure; imagesc(MM'); colorbar;
MM2=MM>0.9;
CC=greedy_cliques(MM2);
figure; imagesc(MM2'); colorbar;

new_templates=zeros(M,T,0);
for j=1:length(CC)
    new_templates=cat(3,new_templates,templates(:,:,CC{j}));
    %new_templates=cat(3,new_templates,zeros(M,T,2));
end;
templates=new_templates;

end

function ret=compute_correlation(v1,v2)
ret=sum(v1.*v2)/sqrt(sum(v1.^2)*sum(v2.^2));
end


function C=greedy_cliques(A)
C={};
while max(A(:))>0
    A0=A; for j=1:size(A,1) A0(j,j)=0; end;
    cc=maximalCliques(A0);
    [maxval,ind0]=max(sum(cc,1));
    if (maxval==1)
        break;
    end;
    ind0=ind0(1);
    C{end+1}=find(cc(:,ind0));
    A(find(cc(:,ind0)),:)=0;
    A(:,find(cc(:,ind0)))=0;
end;
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
