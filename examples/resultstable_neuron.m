% code to generate all results for accuracy table in Neuron paper Chung-Magland et al
% Barnett 7/22/16
%
% datasets on SCDA server have to be at ext_datasets/
%   which links to  /mnt/xfs1/home/ahb/ss_datasets
%   which links to  /mnt/ceph/users/ahb/ss_datasets
%
% parameter (.par) files must also be set up in the data directories.

clear
sorter = @alg_scda_005;
mfile_path = fileparts(mfilename('fullpath'));
ext = [mfile_path,'/../ext_datasets/'];
o_acc.verb = 0;                          % global accuracy options

fout = fopen('resultstable_neuron.txt','w');
fprintf(fout,'dataset  \t M,N \tsamplerate\tCPUtime\tn_gnd\n');   % header
fprintf(fout,'\t\t Qhat extended confusion matrix\n');
fprintf(fout,'\tper neuron...\tk   n_actual   n_good    n_missed    n_falspos \n\n\n');

i=1;   % expt counter (ie dataset)

disp('KAMPFF:')
for neto=1:2          % choose Neto dataset 1 or 2
  d = grab_kampffjuxta_dataset(neto);
  % tmp files will go in /tmp/kampff1output
  o_sort.eleccoordsfile = d.eleccoords;
  o_sort.parfile = [ext,'Kampff/kampff_M32.par'];
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d,o_acc,o_sort);
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  ngnd = 1;  % how many gndtruth neurons
  fprintf(fout,'\nkampff%d:\t %d,%d\t %g\n',neto,d.dims(1),d.dims(2),d.samplerate,cputime{i},ngnd)
  fprintf('\t')
  for j=1:size(Q{i},1), fprintf('%d\t',Q{i}(j,:)), fprintf('\n'); end, fprintf('\n')
  for k=1:ngnd, fprintf('\tk=%d\t%d\t%d\t%d\t%d\n',k,sum(Q{i}(k,:)),Q{i}(k,k),sum(Q{i}(k,:))-Q{i}(k,k),sum(Q{i}(:,k))-Q{i}(k,k)), end
  i=i+1;
end

disp('MARTINEZ:')
for m=1:5             % loop over datasets
  d = grab_martinez2009_dataset(m);
  % tmp files will go in /tmp/kampff1output
  o_acc.parfile = [ext,'MartinezQuiroga2009_sims/martinez.par'];
  t1=tic;
  [fk{i},Q{i},perm{i},info{i}] = accuracy_anysorter_groundtrutheddata(sorter,d,o_acc,os{i});
  cputime{i} = toc(t1); fprintf('done in %.3g sec\n',toc(t1))
  ngnd = 3;  % how many gndtruth neurons
  fprintf(fout,'\nmartinez%d:\t %d,%d\t %g\n',m,d.dims(1),d.dims(2),d.samplerate,cputime{i},ngnd)
  fprintf('\t')
  for j=1:size(Q{i},1), fprintf('%d\t',Q{i}(j,:)), fprintf('\n'); end, fprintf('\n')
  for k=1:ngnd, fprintf('\tk=%d\t%d\t%d\t%d\t%d\n',k,sum(Q{i}(k,:)),Q{i}(k,k),sum(Q{i}(k,:))-Q{i}(k,k),sum(Q{i}(:,k))-Q{i}(k,k)), end
  i=i+1;
end

fclose(fout);
