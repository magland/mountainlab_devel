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

end

