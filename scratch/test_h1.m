function test_h1

close all;

% clips=readmda('franklab/test_hippocampal_01/output/clips.mda');
% features=ms_event_features(clips,6);
% labels=isosplit(features);
% figure; ms_view_clusters(features,labels);

T=60;

fprintf('reading...\n');
X=readmda('franklab/test_hippocampal_01/output/pre3.mda');
detect=readmda('franklab/test_hippocampal_01/output/detect.mda');
clips=ms_extract_clips(X,detect(2,:),T);

clusters=do_sorting(clips);
templates=zeros(size(clips,1),size(clips,2),0);
for j=1:length(clusters)
    templates=cat(3,templates,clusters{j}.template);
end;

figure; 
ms_view_templates(templates);

clips2=zeros(size(clips,1),size(clips,2),0);
labels=[];
for j=1:length(clusters)
    inds=clusters{j}.inds;
    clips2=cat(3,clips2,clips(:,:,inds));
    labels=[labels,ones(1,length(inds))*j];
end;

fprintf('features...\n');
features2=ms_event_features(clips2,6);

figure; ms_view_clusters(features2,labels);

end

function clusters=do_sorting(clips)

[N,T,NC]=size(clips);
Tcenter=ceil((T+1)/2);

clusters={};

isosplit_opts.isocut_threshold=1.2;

fprintf('features...\n');
features=ms_event_features(clips,6);
fprintf('isosplit...\n');
labels=isosplit(features,isosplit_opts);
%figure; ms_view_clusters(features,labels);

K=max(labels);

fprintf('%d clusters found.\n',K);

aa=250;

for k=1:K
    indices_k=find(labels==k);
    clips0=clips(:,:,indices_k);
    if (size(clips0,3)>aa*2)
        peaks=squeeze(max(abs(clips0(:,Tcenter,:)),[],1));
        [~,sort_inds]=sort(peaks);
        inds_to_use=sort_inds(aa:end);
        clips1=clips0(:,:,inds_to_use);
        clusters2=do_sorting(clips1);
        if (length(clusters2)>1)
            for j=1:length(clusters2)
                CC=clusters2{j};
                CC.inds=indices_k(inds_to_use(CC.inds));
                clusters{end+1}=CC;
            end;
        else
            CC.template=mean(clips1,3);
            CC.inds=indices_k;
            clusters{end+1}=CC;
        end;
    else
        CC.template=mean(clips0,3);
        CC.inds=indices_k;
        clusters{end+1}=CC;
    end;
end;

end
