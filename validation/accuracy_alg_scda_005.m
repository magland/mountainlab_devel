% accuracy-measuring function on default dataset, for jfm alg_scda_005
% Barnett 7/28/16

clear
o.sign = -1;      % sorting opts
o.freq_max=inf;      % no roll-off at high freq
o.detectmeth=3;  % 0,3 or 4. 3 seems slightly better than 4
o.detect_threshold = 3.5;
o.clip_size=60;
o.use_mask_out_artifacts=0; o.use_whitening = 1;

oa.verb=1; %2; oa.usepre=1;  % 2 shows spikespy, 3 shows MV (or use below)
[fk Q perm info] = accuracy_anysorter_groundtrutheddata(@alg_scda_005,[],oa,o);
fprintf('alg_scda_005 sort done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');
fprintf('extended confusion matrix:\n')
for i=1:size(Q,1), fprintf('%d\t',Q(i,:)), fprintf('\n'); end

if 0, d=demo_dataset;
  mv.mode='overview2'; mv.raw=d.timeseries; mv.filt=[d.outdir,'/pre2.mda'];
  mv.samplerate=d.samplerate; mv.firings=[d.outdir,'/firings_permed.mda'];
  mountainview(mv);
end
