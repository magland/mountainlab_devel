function mountainview(output_path)

if nargin<1, mountainview_example1; return; end

exe_path=[fileparts(mfilename('fullpath')),'/mountainview/bin/mountainview'];
if (~exist(exe_path,'file'))
    fprintf('Unable to find mountainview executable: %s\n',exe_path);
    fprintf('You need to compile mountainview using Qt4 or Qt5\n');
    fprintf('See readme files for instructions.\n');
    return;
end;

cmd=sprintf('%s --output_path=%s &',exe_path,output_path);
system(cmd);

end

