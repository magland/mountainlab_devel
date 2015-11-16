function mscmd_bandpass_filter(input_path,output_path,opts)

cmd=sprintf('%s bandpass_filter --input=%s --output=%s ',mscmd_exe,input_path,output_path);
cmd=[cmd,sprintf('--samplefreq=%g --freq_min=%g --freq_max=%g',opts.samplefreq,opts.freq_min,opts.freq_max)];

fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end