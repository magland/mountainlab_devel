function ms_setup
mfile_path=fileparts(mfilename('fullpath'));
%addpath([mfile_path,'/mountainlab/matlab']);
%ms_setup_path;
disp('Be sure to run mountainlab_setup.m');

addpath([mfile_path,'/sorting_algs']);
addpath([mfile_path,'/unit_tests']);
addpath([mfile_path,'/validation']);
addpath([mfile_path,'/examples']);
addpath([mfile_path,'/OBSOLETE_mountainview/src/spikespy/matlab']);

end
