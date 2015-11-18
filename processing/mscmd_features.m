function mscmd_features(input_path,detect_path,output_path,opts)

if (nargin<4) opts=struct; end;

cmd=sprintf('%s features --input=%s --detect=%s --output=%s ',mscmd_exe,input_path,detect_path,output_path);
cmd=[cmd,sprintf('--num_features=%d --clip_size=%d --force=1 ',opts.num_features,opts.clip_size)];

fprintf('\n*** FEATURES ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end