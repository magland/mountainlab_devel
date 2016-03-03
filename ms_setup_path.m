function ms_setup_path
mfile_path=fileparts(mfilename('fullpath'));

addpath([mfile_path,'/isosplit']);
addpath([mfile_path,'/mountainview/src/spikespy/matlab']);
addpath([mfile_path,'/msutils']);
addpath([mfile_path,'/processing']);
addpath([mfile_path,'/view']);
addpath([mfile_path,'/unit_tests']);
addpath([mfile_path,'/examples']);
addpath([mfile_path,'/validation']);
end
