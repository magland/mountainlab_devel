function ret=mscmd_exe
exe_fname=sprintf('%s/../bin/mountainsort',fileparts(mfilename('fullpath')));
ret=sprintf('%s %s','LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib',exe_fname);
if (~exist(exe_fname,'file'))
    error('File does not exist: %s\n',exe_fname);
end;
end