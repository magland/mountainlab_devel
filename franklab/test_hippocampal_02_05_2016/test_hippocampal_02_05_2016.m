function test_hippocampal_02_05_2016

close all; drawnow;

mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode'];

tetrode_num=1;
path0=[mfile_path,sprintf('/output_tetrode%d',tetrode_num)];
if ~exist(path0,'dir') mkdir(path0); end;
extract_raw_data(raw_path,path0,tetrode_num);

plausibility_threshold=0.6;
merge_threshold=0.8;
tt_list=4:0.5:15;
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

fprintf('Reading...\n');
clips=readmda([path0,'/clips.mda']);
[M,T,NC]=size(clips);

fprintf('Proj cluster...\n');
label_paths=proj_cluster(clips,tt_list,num_features);

NN=enumerate_paths(label_paths);
print_path_node(NN,'');

% figure;
% aa=label_paths; aa(find(aa(:)==0))=inf;
% for k=1:min(10000,NC)
%     tmp=aa(:,k);
%     plot(1:length(tt_list),tmp+randn(size(tmp))*0.05,'color',rand(1,3));
%     hold on;
% end;


end

function print_path_node(NN,prefix)
if (NN.count<100) return; end;
fprintf('%s%s (%d)\n',prefix,NN.path,NN.count);
while 1
    if (length(NN.children)>1)
        for j=1:length(NN.children)
            print_path_node(NN.children{j},[prefix,'    ']);
        end;
        break
    else
        if (length(NN.children)==1)
            NN=NN.children{1};
        else
            break;
        end;
    end;
end;
end

function NN=enumerate_paths(label_paths,path0)
if nargin<2, path0=[]; end;

NN.count=size(label_paths,2);
NN.path=path0;
NN.children={};
KK=max(label_paths(1,:));
for k=1:KK
    inds_k=find(label_paths(1,:)==k);
    if (length(inds_k)>0)
        if (size(label_paths,1)>1)
            NN_child=enumerate_paths(label_paths(2:end,inds_k),[path0,sprintf('%d',k)]);
        else
            NN_child.path=[path0,sprintf('%d',k)];
            NN_child.count=length(inds_k);
            NN_child.children={};
        end;
        NN.children{end+1}=NN_child;
    end;
end;

end

function label_paths=proj_cluster(clips,tt_list,num_features)

[M,T,NC]=size(clips);

fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

nn=length(tt_list);
label_paths=zeros(nn,NC);

for ii=1:nn
    tt=tt_list(ii);
    inds_tt=find(clip_peaks>=tt);
    fprintf('tt=%d, %d events... ',tt,length(inds_tt));
    clips_tt=clips(:,:,inds_tt);
    fprintf('features... ');
    [FF_tt,subspace_tt]=ms_event_features(clips_tt,num_features);
    fprintf('isosplit... ');
    labels_tt=isosplit(FF_tt);
    K=max(labels_tt);
    fprintf('K=%d\n',K);
    label_paths(ii,inds_tt)=labels_tt;
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
