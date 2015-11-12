function ms_example1_view

addpath([fileparts(mfilename('fullpath')),'/..']);
addpath([fileparts(mfilename('fullpath')),'/../view']);

output_dir_path=[fileparts(mfilename('fullpath')),'/../example1_output'];

mountainview(output_dir_path);
view_cross_correlograms([output_dir_path,'/cross-correlograms.mda'],0);

end
