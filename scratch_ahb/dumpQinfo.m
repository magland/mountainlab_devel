function Qinfo = dumpQinfo(fout,Q,name,d,time,gnd)
% DUMPQINFO   output table format for one dataset validation results against gnd truth
%
% Qinfo = dumpQinfo(fout,Q,name,d,time,gnd)
%
% See: resultstable_neuron.m which uses it
% Barnett 7/25/16, 8/4/16

if nargin==1         % make the header
  fprintf(fout,'dataset  \t M,N \tsamplerate\tCPUtime\tn_gnd\n\n');   % header
  fprintf(fout,'\tQhat extended confusion matrix\n');
  fprintf(fout,'\tk   n_true  frac_missed est_falspos n_missed n_fals\n\n\n');

else                 % output results...
  fprintf(fout,'\n%s:\t %d,%d\t %d\t%.3g\t%d\n\n',name,d.dims(1),d.dims(2),d.samplerate,time,numel(gnd));
  for j=1:size(Q,1), fprintf(fout,'\t%d',Q(j,:)); fprintf(fout,'\n'); end, fprintf(fout,'\n');
  if size(Q,2)<max(gnd)+1            % problem handling
    fprintf('not enough sorted neurons found to match all the gnd-truth ones! size(Q)=[%d,%d]\n',size(Q,1),size(Q,2))
    gnd = gnd(gnd<=size(Q,2)-1);    % just requested that lie within # neurons found
  end
  ntrue = 0; ngood = 0; nmiss = 0; nfals = 0;
  for k=gnd(:)'
    truek = sum(Q(k,:));
    goodk = Q(k,k);
    missk = sum(Q(k,:))-Q(k,k);
    falsk = sum(Q(:,k))-Q(k,k);
    fprintf(fout,'\tk=%d\t%d\t%.3g\t%.3g\t%d\t%d\n',k,truek,missk/truek,falsk/(goodk+falsk),missk,falsk);
    ntrue=ntrue+truek; ngood=ngood+goodk; nmiss=nmiss+missk; nfals=nfals+falsk;
  end
  fprintf(fout,'\ttotal:\t%d\t%.3g\t%.3g\t%d\t%d\n\n',ntrue,nmiss/ntrue,nfals/(ngood+nfals),nmiss,nfals);
  Qinfo.ntrue = ntrue; Qinfo.ngood=ngood; Qinfo.nfals=nfals; Qinfo.nmiss=nmiss;
end
