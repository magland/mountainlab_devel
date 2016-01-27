function mscmd_create_clips_file(input_path,cluster_path,output_path,index_out_path,opts)

if (nargin<5) opts=struct; end;

cmd=sprintf('%s create_clips_file --input=%s --cluster=%s --output=%s --index_out=%s --clip_size=%d ',mscmd_exe,input_path,cluster_path,output_path,index_out_path,opts.clip_size);

fprintf('\n*** CREATE CLIPS FILE ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end
