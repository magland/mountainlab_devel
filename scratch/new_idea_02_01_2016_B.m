function new_idea_02_01_2016_B

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../franklab/test_hippocampal_01_28_2016/tetrode2_output'];

o_detect.threshold=3;
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=10;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=150;

mscmd_detect([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],[mfile_path,'/tmp_clips.mda'],o_extract_clips);

fprintf('Reading...\n');
detect=readmda([mfile_path,'/tmp_detect.mda']);
clips=readmda([mfile_path,'/tmp_clips.mda']);
pre2=readmda([path0,'/pre2.mda']);
[M,T,NC]=size(clips);

%[clips,inds]=random_orthant(clips,6);
%[M,T,NC]=size(clips);
%detect=detect(:,inds);

clips_norms=sqrt(squeeze(sum(sum(clips.^2,1),2)))';
clips_peaks=squeeze(max(clips(:,T/2+1,:),[],1))';

%[~,sort_inds]=sort(clips_peaks);
%clips=clips(:,:,sort_inds);
rand_nums=rand(1,NC);
inds2=find(rand_nums<( (clips_peaks-3)/3 ).^4);
clips=clips(:,:,inds2);
[M,T,NC]=size(clips);

clips=clips./repmat(reshape(clips_norms(inds2),1,1,NC),M,T,1);

fprintf('Features...\n');
FF=ms_event_features(clips,6);
fprintf('isosplit...\n');
labels=isosplit(FF);
figure; ms_view_clusters(FF,labels);
figure; ms_view_templates_from_clips(clips,labels);
%figure; ms_view_clip_clouds(clips,labels);
return;


pow=1;
%clips_pow=clips.*abs(clips).^(pow-1);
clips_pow=clips./repmat(reshape(clips_norms,1,1,NC),M,T,1);
fprintf('Features...\n');
FF=ms_event_features(clips_pow,6);
fprintf('isosplit...\n');
labels=isosplit(FF);
figure; ms_view_clusters(FF,labels);
figure; ms_view_templates_from_clips(clips,labels);
figure; ms_view_templates_from_clips(clips_pow,labels);
figure; ms_view_clip_clouds(clips,labels);

times=detect(2,:);
%spikespy({pre2,times,labels});

clusters=zeros(3,NC);
clusters(1:2,:)=detect;
clusters(3,:)=labels;

writemda(clusters,[mfile_path,'/tmp_clusters.mda']);

writemda(ms_templates(clips,labels),[mfile_path,'/tmp_templates.mda']);

mscmd_cross_correlograms([mfile_path,'/tmp_clusters.mda'],[mfile_path,'/tmp_cross_correlograms.mda'],10000);

[clips_out,clips_index]=ms_create_clips_index(clips,labels);
writemda(clips_out,[mfile_path,'/tmp_clips.mda']);
writemda(clips_index,[mfile_path,'/tmp_clips_index.mda']);

view_params.raw=[path0,'/pre2.mda'];
view_params.clusters=[mfile_path,'/tmp_clusters.mda'];
view_params.templates=[mfile_path,'/tmp_templates.mda'];
view_params.clips=[mfile_path,'/tmp_clips.mda'];
view_params.clips_index=[mfile_path,'/tmp_clips_index.mda'];
view_params.cross_correlograms=[mfile_path,'/tmp_cross_correlograms.mda'];
ms_mountainview(view_params);


return;

clips_norms=sqrt(squeeze(sum(sum(clips.^2,1),2)))';
clips_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clips_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clips_peaks=clips_peaks_pos.*(abs(clips_peaks_pos)>abs(clips_peaks_neg))+clips_peaks_neg.*(abs(clips_peaks_pos)<abs(clips_peaks_neg));

tmp_inds=find(abs(clips_peaks)>5);
tmp_ind=randi(length(tmp_inds));
ind0=tmp_inds(tmp_ind);
clip0=clips(:,:,ind0);
figure; ms_view_templates(clip0);

[closest_clips,resid_clips,closest_inds]=find_closest_clips(clips,clip0,NC);
figure; ms_view_templates(cat(1,closest_clips(:,:,[1:20]),resid_clips(:,:,[1:20])));

pow=4;
aa=closest_clips.*abs(closest_clips).^(pow-1);
bb=clip0.*abs(clip0).^(pow-1);
diffs=aa-repmat(bb,1,1,size(closest_clips,3));
resid_norms=squeeze(sqrt(sum(sum(diffs.^2,1),2)));

randinds=randsample(length(closest_inds),6000);
figure; plot(2:length(resid_norms),resid_norms(2:end),'r');
figure; semilogy(clips_norms(closest_inds(randinds)),resid_norms(randinds),'b.');

return;

scores=compute_scores(clips,clip0,1);

figure; plot(clips_peaks,scores,'b.'); hold on;
plot(clips_peaks(ind0),scores(ind0),'r.','markersize',15);
scores(ind0)=0;

AA=cat(1,clips_peaks,scores);
labels=isosplit(AA);
ms_view_clusters(AA,labels);

end

function [closest_clips,resid_clips,closest_inds]=find_closest_clips(clips,clip0,num)

if nargin<2, clip0=clips; end;

pow=2;
clips_pow=clips.*abs(clips).^(pow-1);
clip0_pow=clip0.*abs(clip0).^(pow-1);

[~,~,NC]=size(clips);
ips=squeeze(sum(sum(clips_pow.*repmat(clip0_pow,1,1,NC),1),2));

clips_pow_norms=sqrt(squeeze(sum(sum(clips_pow.^2,1),2)));
clip0_ref_pow_norm=sqrt(sum(sum(clip0_pow.^2,1),2));
scores1=sqrt(clips_pow_norms.^2+clip0_ref_pow_norm^2-2*ips);

[~,inds_sorted]=sort(scores1,'ascend');
closest_clips=clips(:,:,inds_sorted(1:num));
resid_clips=repmat(clip0,1,1,num)-closest_clips;

closest_inds=inds_sorted(1:num);

end

function [scores,clips_norms]=compute_scores(clips,clips_ref,num)
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
%LARGER_NORMS=max(NORMS1,NORMS2);
%scores1=ips./LARGER_NORMS;
scores1=ips./(NORMS1.*NORMS2);
scores1=ips;
scores1_sorted=sort(scores1,1,'descend');
scores=scores1_sorted(num,:);
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

