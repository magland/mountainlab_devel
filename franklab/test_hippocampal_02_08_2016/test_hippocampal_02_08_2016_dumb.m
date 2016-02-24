function test_hippocampal_02_08_2016


mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode'];

tetrode_num=1;
path0=[mfile_path,sprintf('/output_tetrode%d',tetrode_num)];
if ~exist(path0,'dir') mkdir(path0); end;
extract_raw_data(raw_path,path0,tetrode_num);

plausibility_threshold=0.8;
merge_threshold=0.95;
tt_range=[5,20];
num_tt_steps=20;
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
o_detect.threshold=tt_range(1);
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

fprintf('Reading clips...\n');
clips=readmda([path0,'/clips.mda']);
[M,T,NC]=size(clips);

cluster_test(clips,tt_range,num_tt_steps,num_features);

end

function cluster_test(clips,tt_range,num_tt_steps,num_features)

[M,T,NC]=size(clips);

fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

sorted_peaks=sort(clip_peaks);
inds00=find((sorted_peaks>=tt_range(1))&(sorted_peaks<=tt_range(2)));
sorted_peaks=sorted_peaks(inds00);
tt_list=zeros(1,num_tt_steps);
for jj=1:num_tt_steps
    ind1=max(1,ceil(length(sorted_peaks)*(jj-1)/num_tt_steps+1));
    tt_list(jj)=sorted_peaks(ind1);
end;

clusterings={};
for ii=1:length(tt_list)
    tt0=tt_list(ii); 
    fprintf('tt=%g... ',tt0);
    inds_tt=find(clip_peaks>=tt0);
    CC.tt=tt0;
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
        CC.K=K;
        CC.features=FF_tt;
        CC.clips=clips_tt;
        CC.medoids=zeros(M,T,K);
        for k=1:K
            inds_tt_k=find(labels_tt==k);
            FF_tt_k=FF_tt(:,inds_tt_k);
            [FF_tt_k_medoid,medoid_index]=compute_medoid(FF_tt_k);
            template0=zeros(M,T);
            for jjj=1:size(subspace_tt,3)
                template0=template0+FF_tt_k_medoid(jjj)*subspace_tt(:,:,jjj);
            end;
            CC.medoids(:,:,k)=template0;
            %CC.medoids(:,:,k)=clips_tt(:,:,inds_tt_k(medoid_index));
        end;
    end;
    clusterings{end+1}=CC;
end

clusters={};
for jj=1:length(clusterings)
    CC=clusterings{jj};
    for kk=1:CC.K
        CL.inds=CC.inds(find(CC.labels==kk));
        CL.medoid=CC.medoids(:,:,k);
        CL.tt=CC.tt;
        CL.level=jj;
        clusters{end+1}=CL;
    end;
end;

templates=zeros(M,T,0);
for ii=1:length(clusterings)
    templates=cat(3,templates,clusterings{ii}.medoids);
    templates=cat(3,templates,zeros(M,T,2));
end;


KK=length(clusters);
MM=zeros(KK,KK);
for k1=1:KK
    for k2=1:KK
        val=length(intersect(clusters{k1}.inds,clusters{k2}.inds))/length(clusters{k1}.inds);
        MM(k1,k2)=val;
    end;
end;

parents=zeros(1,KK);
for jj=1:length(clusters)
    inds=find((MM(jj,:)>=0.95)&(1:KK<jj));
    if (length(inds)>0)
        bestind=0;
        bestval=inf;
        for ind=inds
            val=length(clusters{ind}.inds);
            if (val<bestval) bestval=val; bestind=ind; end;
        end;
        parents(jj)=bestind;
    else
        parents(jj)=0;
    end;
end;
figure; treeplot(parents);
[xxx,yyy]=treelayout(parents);
figure; 
for j=1:length(xxx)
    text(xxx(j),yyy(j),sprintf('%d',j));
end;

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

function [m,ind] = compute_medoid(X)
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
