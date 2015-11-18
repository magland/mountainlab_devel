function ret=mscmd_exe
ret=sprintf('%s/../bin/mountainsort',fileparts(mfilename('fullpath')));
if (~exist(ret,'file'))
    error('File does not exist: %s\n',ret);
end;
end