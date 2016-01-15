function test_002

%close all;

ch=1:5;

mfile_path=fileparts(mfilename('fullpath'));
X=readmda(sprintf('%s/../example_data/filt2_white.mda',mfile_path));
X=X(ch,1:1e7);
[M,N]=size(X);

if 1
    opts.detect_interval=200;
    opts.detect_threshold=6;
    opts.clip_size=100;
    [Tpos,Tneg]=ms_detect(X,opts);
    times=sort([Tpos,Tneg]);
    clips=ms_extract_clips(X,times,100);
    writemda(clips,'scratch/clips.mda');
end;
clips=readmda('scratch/clips.mda');
ss_view_waveforms(clips(:,:,1:20));

[~,T,NC]=size(clips);

P=NC;
if 0
    scores=zeros(P,NC);
    norms0=squeeze(sum(sum(clips.^2,1),2))
    for pp=1:P
        fprintf('%d/%d\n',pp,P);
        W0=clips(:,:,pp);
        tmp=clips-repmat(W0,1,1,NC);
        scores(pp,:)=norms0-squeeze(sum(sum(tmp.^2,1),2));
        %scores(scores<0)=0;
    end;
    writemda(scores,'scores.mda');
end;
scores=readmda('scores.mda');
%scores(scores<0)=0;

tmp=zeros(1,NC);
for jj=1:NC
    fprintf('%d/%d\n',jj,NC);
    [~,inds]=sort(scores(:,jj));
    tmp(jj)=scores(inds(end-5),jj);
end;
figure; hist(tmp,1000);
iii=find(tmp>2000);
ss_view_waveforms(clips(:,:,iii(1:30)));
FF=ms_event_features(clips(:,:,iii),3);
labels=isosplit(FF);
ss_view_clusters(FF(1:3,:),labels);

return;

if 1
    best_inds=[];
    for k=1:1000
        disp(k);
        sum_scores=sum(scores,2);
        [~,ind]=max(sum_scores);
        best_inds=[best_inds,ind];
        scores=max(scores,repmat(scores(ind,:),P,1));
    end;
    writemda(best_inds,'best_inds.mda');
end;
best_inds=readmda('best_inds.mda');

ss_view_waveforms(clips(:,:,best_inds(1:30)));

FF=ms_event_features(clips(:,:,best_inds(1:30)),3);
labels=isosplit(FF);
ss_view_clusters(FF(1:3,:),labels);


end
