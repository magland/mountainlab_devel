function test_innermost_shell

clips=readmda('testdata/innermost_shell_clips.mda');
FF1=ms_event_features(clips,12);
labels1=isosplit2(FF1);
figure; ms_view_clusters(FF1,labels1);
title('With whitening at each comparison');

labels2=isosplit2(FF1,struct('whiten_at_each_comparison',0));
figure; ms_view_clusters(FF1,labels2);
title('Without whitening at each comparison');

templates1=ms_templates(clips,labels1);
figure; ms_view_templates(templates1);

K=max(labels1);
for k=1:K
    spikespy({clips(:,:,find(labels1==k)),sprintf('k=%d',k)});
end;

end