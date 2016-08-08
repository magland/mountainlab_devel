% test all detection methods on demodata, compare to known times.
% Barnett 8/4/16

clear
d = demo_dataset;
f = readmda(d.truefirings);
Y = readmda(d.timeseries);
%Y = ms_normalize_channels(Y);
Y = ms_whiten(Y);
od=[]; od.clip_size = 60; od.detect_interval=10; od.sign=-1;
od.detect_threshold = 3.5;
od.individual_channels = 1;   % 1 actually detects better than 0

meths= [0   1   3 3  4 4];
betas= [nan nan 1 10 1 10];   % upsampling if relevant
for i=1:numel(meths)
  meth=meths(i); od.beta = betas(i); fprintf('\nmeth=%d,beta=%d..........\n',meth,od.beta)
  if meth==0
    outdir = tempdir; detectfile = [outdir,'/detect.mda'];  % ie run_mscmd_detect
    mscmd_detect(pathify32(Y),detectfile,od);
    firings = readmda(detectfile);
    times = firings(2,:);
  elseif meth==1, times = ms_detect(Y,od);       % obsolete, never individ_chan
  elseif meth==3, [times chans] = ms_detect3(Y,od);
  elseif meth==4, [times chans] = ms_detect4(Y,od);
  end
  fprintf('\t %d events detected\n',numel(times))
  Ttrue = f(2,:); Ltrue = f(3,:);
  oo.verb = 0;
  [perm Q acc t]=times_labels_accuracy(Ttrue,Ltrue,times,1+0*times,oo);
end
%Q
%spikespy({Y,Ttrue,Ltrue,'true'},{Y,times,1+0*times,'det'},{Y,t.tmiss,1+0*t.tmiss,'tmiss'});
%display_confusion(Q,'true','just detection');
