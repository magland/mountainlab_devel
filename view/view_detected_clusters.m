function view_detected_clusters(ch,opts,sign)

posneg='pos';
if (sign<0) posneg='neg'; end;
fname=sprintf('%s%s_%d.mda',opts.cluster_features_prefix,posneg,ch);
try
    FF=readmda(fname);
catch
    error(sprintf('Unable to read file %s. Perhaps you need to run the processing.',fname));
end
fname=sprintf('%s%s_%d.mda',opts.cluster_labels_prefix,posneg,ch);
labels=readmda(fname);
opts.create_figure=0;
ss_view_clusters(FF,labels,opts);

end