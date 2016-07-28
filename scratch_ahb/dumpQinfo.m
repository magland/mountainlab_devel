function dumpQinfo(fout,Q,name,d,time,gnd)
% DUMPQINFO   output table format for one dataset validation results against gnd truth
%
% dumpQinfo(fout,Q,name,d,time,gnd)
%
% See: resultstable_neuron.m which uses it
% Barnett 7/25/16

fprintf(fout,'\n%s:\t %d,%d\t %d\t%.3g\n\n',name,d.dims(1),d.dims(2),d.samplerate,time);
for j=1:size(Q,1), fprintf(fout,'\t%d',Q(j,:)); fprintf(fout,'\n'); end, fprintf(fout,'\n');
if size(Q,2)<max(gnd)+1            % problem handling
  fprintf('not enough sorted neurons found to match all the gnd-truth ones! size(Q)=[%d,%d]\n',size(Q,1),size(Q,2))
  gnd = gnd(gnd<=size(Q,2)-1);    % just requested that lie within # neurons found
end
nacc = 0; ngood = 0; nmiss = 0; nfals = 0;
for k=gnd(:)'
  fprintf(fout,'\tk=%d\t%d\t%d\t%d\t%d\n',k,sum(Q(k,:)),Q(k,k),sum(Q(k,:))-Q(k,k),sum(Q(:,k))-Q(k,k));
  nacc=nacc+sum(Q(k,:)); ngood=ngood+Q(k,k); nmiss=nmiss+sum(Q(k,:))-Q(k,k);
  nfals=nfals+sum(Q(:,k))-Q(k,k);
end
fprintf(fout,'\ttotal:\t%d\t%d\t%d\t%d\n\n',nacc,ngood,nmiss,nfals);
