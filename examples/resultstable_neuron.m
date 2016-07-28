% code to generate all results for accuracy table in Neuron paper Chung-Magland et al
% Barnett 7/22/16
%
% datasets on SCDA server have to be at ext_datasets/
%   which links to  /mnt/xfs1/home/ahb/ss_datasets
%   which links to  /mnt/ceph/users/ahb/ss_datasets
%
% parameter (.par) files must also be set up in the data directories, pointed to here.
% Needs: dumpQinfo.m

clear
sorter = @alg_scda_005_js;                    % wrapper to js sorter chain
mfile_path = fileparts(mfilename('fullpath'));
ext = [mfile_path,'/../ext_datasets/'];
o_acc.verb = 0;                          % global accuracy options

fout = fopen('resultstable_neuron.txt','w');
fprintf(fout,'dataset  \t M,N \tsamplerate\tCPUtime\tn_gnd\n\n');   % header
fprintf(fout,'\tQhat extended confusion matrix\n');
fprintf(fout,'\tk   n_actual   n_good    n_missed    n_falspos \n\n\n');

i=1;   % expt counter (ie dataset)

if 0
disp('KAMPFF:')
% Compare to text p.10 of Neto '16 biorxiv preprint.
for neto=1:2          % choose Neto dataset 1 or 2  (1 gives only K=12 good, 2 K=24ish)
  d{i} = grab_kampffjuxta_dataset(neto);
  % tmp files will go in /tmp/kampff1output
  o_sort.eleccoordsfile = d{i}.eleccoords;
  o_sort.parfile = [ext,'Kampff/kampff_M32.par'];
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d{i},o_acc,o_sort);
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  dumpQinfo(fout,Q{i},sprintf('kampff%d',neto),d{i},cputime{i},1);  % 1 gnd-truth neuron
  i=i+1;
end
end


if 0
disp('MARTINEZ:')
% compare to: ~/spikesorting/martinezresults.txt which was for jfm_april_sort
% and compare tables in Martinez '09 paper.
% Params note: sign=+1 and need switch off detectability thresh
for m=1:5             % loop over datasets
  d{i} = grab_martinez2009_dataset(m);    % tmp files will go in /tmp/martinez1output etc
  o_sort = [];       % remove the eleccoords
  o_sort.parfile = [ext,'MartinezQuiroga2009_sims/martinez.par'];
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d{i},o_acc,o_sort);
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  dumpQinfo(fout,Q{i},sprintf('martinez%d',m),d{i},cputime{i},2:3);  % 2/3 gnd-truth neurons are the single-unit (k=1 is MUA)
  i=i+1;
end
end


if 1
  disp('HARRIS2000:')
% and compare tables in Chaitu etc.
% Params note: sign=-1 and need switch off detectability thresh
  d{i} = grab_harris2000_dataset;    % tmp files will go in /tmp/harrisoutput
  o_sort = [];       % remove the eleccoords
  o_sort.parfile = [ext,'Harris2000/d5331/harris2000_5331.par'];
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d{i},o_acc,o_sort);
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  dumpQinfo(fout,Q{i},'harris2000',d{i},cputime{i},1);
  i=i+1;
end











fclose(fout);

% Generic sorting view of the i'th thing that was run:
i=1; mv.raw=d{i}.timeseries; mv.filt=[d{i}.outdir '/pre2.mda']; mv.firings=[d{i}.outdir '/firings.mda']; mv.samplerate=d{i}.samplerate; mountainview(mv);
