function mscmd_detect(input_path,output_path,opts)

if (nargin<3) opts=struct; end;

cmd=sprintf('%s detect --input=%s --output=%s ',mscmd_exe,input_path,output_path);
cmd=[cmd,sprintf('--inner_window_width=%d --outer_window_width=%d ',opts.inner_window_width,opts.outer_window_width)];
cmd=[cmd,sprintf('--threshold=%g ',opts.threshold)];

fprintf('\n*** DETECT ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end