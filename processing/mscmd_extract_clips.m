function mscmd_extract_clips(input_path,cluster_path,output_path,index_out_path,opts)

if (nargin<5) opts=struct; end;

cmd=sprintf('%s extract_clips --input=%s --cluster=%s --output=%s --index_out=%s --clip_size=%d ',mscmd_exe,input_path,cluster_path,output_path,index_out_path,opts.clip_size);

fprintf('\n*** EXTRACT CLIPS ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end
