function mscmd_extract_clips(input_path,detect_path,output_clips_path,opts)

if (nargin<4) opts=struct; end;

cmd=sprintf('%s extract_clips --input=%s --detect=%s --output_clips=%s --clip_size=%d ',mscmd_exe,input_path,detect_path,output_clips_path,opts.clip_size);

fprintf('\n*** EXTRACT CLIPS ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end
