function mscmd_normalize_channels(input_path,output_path)

cmd=sprintf('%s normalize_channels --input=%s --output=%s ',mscmd_exe,input_path,output_path);

fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end