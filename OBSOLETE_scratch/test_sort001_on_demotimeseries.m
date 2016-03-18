function test_sort001_on_demotimeseries

close all;

%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../unit_tests/demo_data/demotimeseries.mda'];
if ~exist(raw_path), writedemotimeseries; end

mfile_path=fileparts(mfilename('fullpath'));
output_path=[mfile_path,'/output_sort001_on_demotimeseries'];
if ~exist(output_path,'dir') mkdir(output_path); end;

%%%% Sort
sort_opts.shell.merge_threshold=0.5;
sort_opts.shell.section_increment=1;
sort_opts.shell.num_sections_per_shell=2;
sort_opts.shell.min_section_count=200;
sort_opts.detect.detect_interval=50;
sort_opts.detectibility_threshold=3;
sort_opts.isosplit.isocut_threshold=1.5;
sort_opts.isosplit.K_init=50;
[firings_path,pre_path]=sort_001(raw_path,output_path,sort_opts);

%%%% View output
mv.mode='overview2';
mv.raw=pre_path;
mv.firings=firings_path;
ms_mountainview(mv);


end