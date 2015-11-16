function mscmd_extract(input_path,output_path,opts)

if (~isfield(opts,'threshold')) opts.threshold=0; end;

cmd=sprintf('%s extract --input=%s --output=%s ',mscmd_exe,input_path,output_path);

channels_str='';
for ii=1:length(opts.channels)
    if (ii>1) channels_str=[channels_str,',']; end;
    channels_str=[channels_str,sprintf('%d',opts.channels(ii))];
end;

cmd=[cmd,sprintf('--num_channels=%d --channels=%s --t1=%d --t2=%d --threshold=%g',opts.num_channels,channels_str,opts.t1,opts.t2,opts.threshold)];

fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end
