% code to generate all results for accuracy table in Neuron paper Chung-Magland et al
% Barnett 7/22/16
%
% datasets on SCDA server have to be at ext_datasets/
%   which links to  /mnt/xfs1/home/ahb/ss_datasets
%   which links to  /mnt/ceph/users/ahb/ss_datasets
%
% parameter (.par) files must also be set up in the data directories, pointed to here.
% Needs: dumpQinfo.m
% Barnett late July - 8/4/16

clear; addpath scratch_ahb
%sorter = @alg_scda_005_js;                 % wrapper to js sorter chain
sorter = @alg_scda_005;                    % current matlab sorter
go.detectmeth = 0;    % 0 (3 no better).  % general sorting opts
mfile_path = fileparts(mfilename('fullpath'));
ext = [mfile_path,'/../ext_datasets/'];
o_acc.verb = 0;                          % global accuracy options

fout = fopen('~/spikesorting/franklab_paper/resultstable_neuron_harr.txt','w');
dumpQinfo(fout);   % does the header
i=1;   % expt counter (ie dataset)

if 1
  disp('HARRIS2000:')
% and compare tables in Chaitu etc.
% Params note: sign=-1 and need switch off detectability thresh
  d{i} = grab_harris2000_dataset;    % tmp files will go in /tmp/harrisoutput
  os = go;
  %os.parfile = [ext,'Harris2000/d5331/harris2000_5331.par']; % for js
  os.sign=-1, os.clip_size=20; os.detect_threshold=3.5; os.freq_min=300;
  os.freq_max=inf; os.use_whitening=0; os.use_mask_out_artifacts=0;
  os.num_fea = 20;
  %os.detectmeth = 0;  % 3,4 no better
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d{i},o_acc,os);
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  Qi{i} = dumpQinfo(fout,Q{i},'harris2000',d{i},cputime{i},1);
  i=i+1;
end


if 0
disp('MARTINEZ:')
% compare to: ~/spikesorting/martinezresults.txt which was for jfm_april_sort
% and compare tables in Martinez '09 paper.
% Params note: sign=+1 and need switch off detectability thresh. M=1
for m=1:5             % loop over datasets
  d{i} = grab_martinez2009_dataset(m);    % tmp files will go in /tmp/martinez1output etc
  os = go;
  %os.parfile = [ext,'MartinezQuiroga2009_sims/martinez.par'];  % for js
  os.clip_size=40; os.sign=1; os.freq_min=300; os.freq_max=10000;
  os.use_whitening=0; os.use_mask_out_artifacts=0; os.detect_threshold=3.0;
  %os.detectmeth = 3;   % *** alarmingly, set#3 wrongly merges 2,3 if det=3 !
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d{i},o_acc,os);
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  Qi{i} = dumpQinfo(fout,Q{i},sprintf('martinez%d',m),d{i},cputime{i},2:3);  % 2/3 gnd-truth neurons are the single-unit (k=1 is MUA)
  i=i+1;
end
end


if 0
disp('NEUROCUBE:')
% and compare tables in Camunas-Mesa '14 paper.
% Params note: sign=+1 and need switch off detectability thresh
  d{i} = grab_neurocube_dataset(2);    % tmp files will go in /tmp/neurocubeoutput
  truefirings = arrayify(d{i}.truefirings); Ktrue = max(truefirings(3,:));
  os = go;
  os.clip_size=50; os.sign=-1; os.freq_min=0; os.freq_max=inf;
  os.use_whitening=0; os.use_mask_out_artifacts=0; os.detect_threshold=3.0;
  os.num_fea = 20;    % default 10; slows it down
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d{i},o_acc,os);
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  
  kk = find(fk{i}>0.5);
  % *** select just the large f_k neurons...
  Qi{i} = dumpQinfo(fout,Q{i},'neurocube',d{i},cputime{i},1:Ktrue);
  
  i=i+1;
end


if 0
disp('KAMPFF/NETO:')
% Compare to text p.10 of Neto '16 biorxiv preprint. 1 JC elec.
for neto=1:2          % choose Neto dataset 1 or 2  (1 gives only K=12 good, 2 K=24ish)
  d{i} = grab_kampffjuxta_dataset(neto);
  % tmp files will go in /tmp/kampff1output
  %os.parfile = [ext,'Kampff/kampff_M32.par']; % for js
  os = go;
  os.eleccoordsfile = d{i}.eleccoords; os.adj_radius=1.5;
  os.clip_size=60; os.sign=-1;   % sign crucial
  os.freq_min=300;os.freq_max=10000;
  os.use_whitening=1; os.use_mask_out_artifacts=0; os.detect_threshold=3.5;
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d{i},o_acc,os);
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  Qi{i} = dumpQinfo(fout,Q{i},sprintf('kampff%d',neto),d{i},cputime{i},1);  % 1 gnd-truth neuron
  i=i+1;
end
end


fclose(fout);

% Generic sorting view of the i'th thing that was run:
%i=1; mv.raw=d{i}.timeseries; mv.filt=[d{i}.outdir '/pre2.mda']; mv.firings=[d{i}.outdir '/firings_permed.mda']; mv.samplerate=d{i}.samplerate; mountainview(mv);

% spikespy comparison of sorted vs truth:
%i=1; Y=readmda(d{i}.timeseries); F=readmda([d{i}.outdir '/firings_permed.mda']); Ft=readmda(d{i}.truefirings); spikespy({Y,F(2,:),F(3,:),'sorted'},{Y,Ft(2,:),Ft(3,:),'truth'});

% gndtruth in MV:
%i=1; mv.raw=d{i}.timeseries; mv.filt=[d{i}.outdir '/pre2.mda']; mv.firings=d{i}.truefirings; mv.samplerate=d{i}.samplerate; mountainview(mv);

save ~/spikesorting/franklab_paper/resultstable_neuron_harr.mat

%i=1; display_confusion(Q{1},'true','sorted');
