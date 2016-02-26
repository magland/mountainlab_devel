function mscmd_remove_duplicates(firings_in_path,firings_out_path,opts)

if (nargin==0) test_mscmd_remove_duplicates; return; end;

if (nargin<3) opts=struct; end;

if (~isfield(opts,'max_dt')) opts.max_dt=6; end;
if (~isfield(opts,'overlap_threshold')) opts.overlap_threshold=0.25; end;

cmd=sprintf('%s remove_duplicates --firings_in=%s --firings_out=%s --max_dt=%d --overlap_threshold=%g',mscmd_exe,firings_in_path,firings_out_path,...
    opts.max_dt,opts.overlap_threshold);

fprintf('\n*** REMOVE_DUPLICATES ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end

function test_mscmd_remove_duplicates

K=6;
num_events=1000;

firings=zeros(3,0);
t=1;
for i=1:num_events
    k=randi(K);
    firings(:,end+1)=[0,t,k];
    if (k==3)
        firings(:,end+1)=[0,t+2,k+1];
    end;
    if (k==5)
        firings(:,end+1)=[0,t-4,k+1];
    end;
    t=t+randi(50);
end;

writemda(firings,'tmp_firings_in.mda');
mscmd_remove_duplicates('tmp_firings_in.mda','tmp_firings_out.mda');
firings2=readmda('tmp_firings_out.mda');
size(firings)
size(firings2)
max(firings2(3,:))

end

