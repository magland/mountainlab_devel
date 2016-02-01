function new_idea_01_30_2016

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0='franklab/test_hippocampal_01_28_2016/tetrode2_output';

o_detect.threshold=3;
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=100;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=100;

mscmd_detect([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],[mfile_path,'/tmp_clips.mda'],o_extract_clips);

clips=readmda([mfile_path,'/tmp_clips.mda']);
[M,T,NC]=size(clips)

use_it=zeros(1,NC);

for j=1:300
    fprintf('j=%d\n',j);
    [clips0,inds0]=random_orthant(clips,6);
    NC0=size(clips0,3);

    [A,norms]=compute_scores_matrix(clips0);
    Asorted=sort(A,1,'descend');
    aa=Asorted(5,:);
    %figure; plot(norms,aa,'b.');
    inds1=find(aa>0.5);
    use_it(inds0(inds1))=use_it(inds0(inds1))+1;
    %figure; imagesc(A); colorbar;
    fprintf('Using %d of %d clips.\n',count(se_it),length(use_it));
end;

fprintf('Using %d of %d clips.\n',sum(use_it),length(use_it));

clips1=clips(:,:,find(use_it));
FF=ms_event_features(clips1,6);
labels=isosplit(FF);
figure; ms_view_clusters(FF,labels);
figure; ms_view_templates_from_clips(clips1,labels);
figure; ms_view_templates(get_example_clips(clips1,labels,5));

return;

Asorted=sort(abs(A),1,'descend');
figure; imagesc(Asorted); colorbar;
aa=Asorted(5,:);
figure; plot(norms,aa,'b.');
inds=find(aa>0.5); 
clips1=clips0(:,:,inds);

FF=ms_event_features(clips1,3);
labels=isosplit(FF);
figure; ms_view_clusters(FF,labels);
figure; ms_view_templates_from_clips(clips1,labels);
figure; ms_view_templates(get_example_clips(clips1,labels,10));

% [evecs,evals]=eig(A); evals=diag(evals);
% figure; plot(evals);
% evec0=evecs(:,end);
% clip0=sum(clips0.*repmat(reshape(evec0,1,1,NC0),M,T,1),3);
% ms_view_templates(clip0);
% figure; hist(evec0,1000);

return;

FF=ms_event_features(clips0,3);
labels=isosplit(FF);
figure; ms_view_clusters(FF,labels);
figure; ms_view_templates_from_clips(clips0,labels);

return;

[M,T,NC]=size(clips0)
clip_norms=reshape(sqrt(sum(sum(clips0.^2,1),2)),1,NC);

numrand=1000000;
num_best=1000;
adj=0;

scores=zeros(1,numrand);
indpairs=randi(NC,2,numrand);
while 1
    tmp1=find(indpairs(1,:)==indpairs(2,:));
    if (length(tmp1)==0) break; end;
    indpairs(:,tmp1)=randi(NC,2,length(tmp1));
end

ips=reshape(sum(sum(clips0(:,:,indpairs(1,:)).*clips0(:,:,indpairs(2,:)),1),2),1,numrand);
larger_norms=max(clip_norms(indpairs(1,:)),clip_norms(indpairs(2,:)));
scores=ips./larger_norms.^2;

%aa=find(larger_norms<40);
%larger_norms=larger_norms(aa);
%scores=scores(aa);
%indpairs=indpairs(:,aa);


% figure;
% subplot(2,1,1);
% plot(larger_norms,scores+larger_norms*adj,'b.');
% subplot(2,1,2);
% plot(larger_norms,scores-larger_norms*adj,'b.');

[~,inds0]=sort(scores,'descend');
best_inds1=indpairs(1,inds0(1:num_best));
best_inds2=indpairs(2,inds0(1:num_best));

inds_to_use=zeros(1,NC);
inds_to_use(best_inds1)=1;
inds_to_use(best_inds2)=1;
best_inds=find(inds_to_use);

best_clips=clips0(:,:,best_inds);
num_best=length(best_inds);
selection=randsample(num_best,30);
figure; ms_view_templates(best_clips(:,:,selection));
FF=ms_event_features(best_clips,3);
labels=isosplit(FF);
figure; ms_view_clusters(FF,labels);

figure; ms_view_templates_from_clips(best_clips,labels);

return;

[bandwidth,density,X,Y]=kde2d(cat(1,larger_norms,scores)',2^8,[0,-1],[max(larger_norms),1]); 

figure; imagesc(density(end:-1:1,:));
%figure; surf(X,Y,density);

density2=log(density)-log(density(end:-1:1,:));
figure; imagesc(max(-3,min(3,density2(end:-1:1,:)))); colorbar;

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

function [A,clip_norms]=compute_scores_matrix(clips)
[M,T,NC]=size(clips);
clip_norms=reshape(sqrt(sum(sum(clips.^2,1),2)),1,NC);
[ii1,ii2]=ndgrid(1:NC,1:NC);
largest_clip_norms=reshape(max(clip_norms(ii1(:)),clip_norms(ii2(:))),NC,NC);
ips=zeros(NC,NC);
for j=1:NC
    if (mod(j,50)==0) disp(j); end;
    CC=clips(:,:,j);
    ips(j,:)=reshape(sum(sum(repmat(CC,1,1,NC).*clips,1),2),1,NC);
end;
A=ips./largest_clip_norms.^2;
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

