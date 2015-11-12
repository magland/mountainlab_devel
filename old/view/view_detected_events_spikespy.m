function view_detected_events_spikespy(ch,opts,data)

fname_times=sprintf('%s%d.mda',opts.detect_times_prefix,ch);
fname_labels=sprintf('%s%d.mda',opts.cluster_labels_prefix,ch);
X=data.X;
try
    T=readmda(fname_times);
    L=readmda(fname_labels);
catch
    error(sprintf('Unable to read file %s or %s. Perhaps you need to run the processing.',fname_times,fname_labels));
end;

AM=readmda(opts.adjacency);
AM(ch,ch)=0;
channel_inds=[ch;find(AM(:,ch))];
X=X(channel_inds,:);

if (length(X(:))>1e7)
    fprintf('Launching spikespy for large array, this could take some time on the first run...\n');
    fprintf('Also note that I am putting large files into to your temporary directory.\n');
end;
spikespy({X,T,L});

end
