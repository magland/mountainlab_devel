function mscmd_detect(input_path,output_path,opts)

if (nargin<3) opts=struct; end;

if ~isfield(opts,'individual_channels') opts.individual_channels=1; end;

cmd=sprintf('%s detect --input=%s --output=%s ',mscmd_exe,input_path,output_path);
cmd=[cmd,sprintf('--inner_window_width=%d --outer_window_width=%d --individual_channels=%d ',opts.inner_window_width,opts.outer_window_width,opts.individual_channels)];
cmd=[cmd,sprintf('--threshold=%g ',opts.threshold)];

fprintf('\n*** DETECT ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end