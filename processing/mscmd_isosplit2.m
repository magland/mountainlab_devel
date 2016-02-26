function mscmd_isosplit2(input_path,out_labels_path,opts)

if (nargin==0) test_mscmd_isosplit2; return; end;

if (nargin<3) opts=struct; end;

input_path=mktmpfile(input_path);

if (~isfield(opts,'isocut_threshold')) opts.isocut_threshold=1.5; end;
if (~isfield(opts,'K_init')) opts.K_init=30; end;

cmd=sprintf('%s isosplit2 --input=%s --labels=%s --isocut_threshold=%g --K_init=%d',mscmd_exe,input_path,out_labels_path,...
    opts.isocut_threshold,opts.K_init);

fprintf('\n*** ISOSPLIT2 ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end

function test_mscmd_isosplit2

M=3;
N=500;
X=randn(M,N);
X(1,300:end)=X(1,300:end)+6;

writemda(X,'tmp_X.mda');
mscmd_isosplit2('tmp_X.mda','tmp_labels.mda');
labels=readmda('tmp_labels.mda');

ms_view_clusters(X,labels);

end