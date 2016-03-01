function ds001_sort_all

mfile_path=fileparts(mfilename('fullpath'));

path0=[mfile_path,'/output_tetrode1'];
[firings_path1,pre_path1]=ds001_sort([path0,'/pre0.mda'],path0);
mv.raw=pre_path1; mv.firings=firings_path1;
ms_mountainview(mv);

path0=[mfile_path,'/output_tetrode2'];
[firings_path1,pre_path1]=ds001_sort([path0,'/pre0.mda'],path0);
mv.raw=pre_path1; mv.firings=firings_path1;
ms_mountainview(mv);

path0=[mfile_path,'/output_ms11d45'];
sort_opts=struct;
sort_opts.adjacency_matrix=[path0,'/adjacency_matrix.mda'];
[firings_path1,pre_path1]=ds001_sort([path0,'/pre0.mda'],path0,sort_opts);
mv.raw=pre_path1; mv.firings=firings_path1;
ms_mountainview(mv);

path0=[mfile_path,'/output_synth_EJ_K7'];
[firings_path1,pre_path1]=ds001_sort([path0,'/pre0.mda'],path0);
mv.raw=pre_path1; mv.firings=firings_path1;
ms_mountainview(mv);


end

