function new_idea_02_01_2016

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../franklab/test_hippocampal_01_28_2016/tetrode2_output'];

o_detect.threshold=4;
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=50;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=150;

mscmd_detect([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],[mfile_path,'/tmp_clips.mda'],o_extract_clips);

fprintf('Reading...\n');
detect=readmda([mfile_path,'/tmp_detect.mda']);
clips=readmda([mfile_path,'/tmp_clips.mda']);
pre2=readmda([path0,'/pre2.mda']);
[M,T,NC]=size(clips);

clips_norms=sqrt(squeeze(sum(sum(clips.^2,1),2)))';
clips_peaks=squeeze(max(clips(:,T/2+1,:).^2,[],1))';

scores=compute_scores_via_partition(clips);
figure; plot(clips_norms,scores,'b.');
xlabel('Norms');
ylabel('Scores');
figure; plot(clips_peaks,scores,'m.');
xlabel('Peaks');
ylabel('Scores');

score_cutoff=5;
peak_cutoff=4;
inds0=find((clips_peaks>peak_cutoff)&(scores>score_cutoff));
fprintf('Using %d of %d clips...\n',length(inds0),NC);
clips0=clips(:,:,inds0);
[~,~,NC0]=size(clips0);

FF0=ms_event_features(clips0,6);
labels0=isosplit(FF0);
K0=max(labels0);

[clips0_out,clips_index]=ms_create_clips_index(clips0,labels0);
writemda(clips0_out,[mfile_path,'/tmp_clips0.mda']);
writemda(clips_index,[mfile_path,'/tmp_clips0_index.mda']);

figure; ms_view_clusters(FF0,labels0);
templates=ms_templates(clips0,labels0);
writemda(templates,[mfile_path,'/tmp_templates.mda']);
figure; ms_view_templates_from_clips(clips0,labels0);

times0=detect(2,inds0);
spikespy({pre2,times0,labels0});

clusters=zeros(3,NC0);
clusters(1:2,:)=detect(:,inds0);
clusters(3,:)=labels0;
writemda(clusters,[mfile_path,'/tmp_clusters.mda']);
mscmd_cross_correlograms([mfile_path,'/tmp_clusters.mda'],[mfile_path,'/tmp_cross_correlograms.mda'],10000);
view_params.raw=[path0,'/pre2.mda'];
view_params.clusters=[mfile_path,'/tmp_clusters.mda'];
view_params.templates=[mfile_path,'/tmp_templates.mda'];
view_params.clips=[mfile_path,'/tmp_clips0.mda'];
view_params.clips_index=[mfile_path,'/tmp_clips0_index.mda'];
view_params.cross_correlograms=[mfile_path,'/tmp_cross_correlograms.mda'];
ms_mountainview(view_params);

% for k=1:K0
%     clips_k=clips0(:,:,find(labels0==k));
%     fprintf('Computing scores for k=%d\n',k);
%     [scores_all_k,norms_all_k]=compute_scores(clips,clips_k);
%     [scores_k,norms_k]=compute_scores(clips_k);
%     figure; plot(norms_all_k,scores_all_k,'k.',norms_k,scores_k,'r.');
%     xlabel('Norms');
%     ylabel('Scores');
%     title(sprintf('Scores for k=%d\n',k));
% end;

for aa=1:1
    example_clips=get_example_clips(clips0,labels0,6);
    figure; ms_view_templates(example_clips);
end;

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
num_splits=8;
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
[~,~,NC]=size(clips);
[~,~,NC_ref]=size(clips_ref);
ips=zeros(NC_ref,NC);
for j=1:NC_ref
    ips(j,:)=squeeze(sum(sum(clips.*repmat(clips_ref(:,:,j),1,1,NC),1),2));
end;
clips_norms=sqrt(squeeze(sum(sum(clips.^2,1),2)));
clips_ref_norms=sqrt(squeeze(sum(sum(clips_ref.^2,1),2)));
[NORMS1,NORMS2]=ndgrid(clips_ref_norms,clips_norms);
LARGER_NORMS=max(NORMS1,NORMS2);
scores1=ips./LARGER_NORMS;
%scores1=ips./(NORMS1.*NORMS2);
scores1_sorted=sort(scores1,1,'descend');
scores=scores1_sorted(20,:);
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



