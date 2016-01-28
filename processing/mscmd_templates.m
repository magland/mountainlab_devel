function mscmd_templates(input_path,cluster_path,output_path,opts)

if (nargin<4) opts=struct; end;

cmd=sprintf('%s templates --input=%s --clusters=%s --output=%s --clip_size=%d ',mscmd_exe,input_path,cluster_path,output_path,opts.clip_size);

fprintf('\n*** TEMPLATES ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end