function view_detected_waveforms(ch,opts)

fname=sprintf('%s%d.mda',opts.cluster_waveforms_prefix,ch);
try
    X=readmda(fname);
catch
    error(sprintf('Unable to read file %s. Perhaps you need to run the processing.',fname));
end
opts.create_figure=0;
ss_view_waveforms(X,opts);

end