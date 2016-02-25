function test_sort_probe_data_02_24_2016

close all;

rng(1);

mfile_path=fileparts(mfilename('fullpath'));
addpath(sprintf('%s/../sort_002_multichannel',mfile_path));
path0=sprintf('%s/output',mfile_path);
if (~exist(path0,'dir')) mkdir(path0); end;
datfile_path=sprintf('%s/../raw/ms11d45.dat',mfile_path);
locations=get_frank_lab_locations;
adjacency_matrix=ms_adjacency_matrix(locations,2);

o_extract.num_channels=72;
o_extract.channels=[37:52,68,69];
o_extract.t1=0; o_extract.t2=19e6;
addpath(sprintf('%s/../../processing/old',mfile_path));
mscmd_extract(datfile_path,[path0,'/raw.mda'],o_extract);

sort_opts.adjacency_matrix=adjacency_matrix;
sort_opts.merge_threshold=0.5;
%sort_opts.debug_channels=[2];
sort_002_multichannel([path0,'/raw.mda'],path0,sort_opts);

%%%% View output
mv.mode='overview2';
mv.raw=[path0,'/pre2.mda'];
mv.firings=[path0,'/firings.mda'];
ms_mountainview(mv);

end

function L=get_frank_lab_locations

L=[...
0.0,0.0;...
-0.5,1.0;...
0.5,1.0;...
-1.0,2.0;...
1.0,2.0;...
-1.0,3.0;...
1.0,3.0;...
-1.0,4.0;...
1.0,4.0;...
-1.0,5.0;...
1.0,5.0;...
-1.0,6.0;...
1.0,6.0;...
-1.0,7.0;...
1.0,7.0;...
-1.0,8.0;...
1.0,8.0;...
-1.0,9.0;...
];

end