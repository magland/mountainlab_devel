function mscmd_outlier_scores_v1(raw_path,firings_in_path,firings_out_path,opts)

if (nargin==0) test_mscmd_outlier_scores_v1; return; end;

if (nargin<4) opts=struct; end;

if (~isfield(opts,'clip_size')) opts.clip_size=100; end;
if (~isfield(opts,'min_shell_size')) opts.min_shell_size=100; end;
if (~isfield(opts,'shell_increment')) opts.shell_increment=0.5; end;

cmd=sprintf('%s outlier_scores_v1 --raw=%s --firings_in=%s --firings_out=%s --clip_size=%d --min_shell_size=%d --shell_increment=%g',mscmd_exe,raw_path,firings_in_path,firings_out_path,...
    opts.clip_size,opts.min_shell_size,opts.shell_increment);

fprintf('\n*** OUTLIER_SCORES_V1 ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end

function test_mscmd_outlier_scores_v1
mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../demo/demo_sort_002/output_tetrode1'];
mscmd_outlier_scores_v1([path0,'/pre2.mda'],[path0,'/firings.mda'],'tmp_firings.mda');
mv.pre=[path0,'/pre2.mda'];
mv.firings='tmp_firings.mda';
mv.sampling_freq=30000;
ms_mountainview(mv);
end
