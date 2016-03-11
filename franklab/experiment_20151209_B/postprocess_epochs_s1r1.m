function postprocess_epochs_s1r1

mfile_path=fileparts(mfilename('fullpath'));
pre2a=readmda([mfile_path,'/output_s1/pre2.mda']);
pre2b=readmda([mfile_path,'/output_r1/pre2.mda']);
firings=readmda([mfile_path,'/output_s1r1/firings.mda']);

times=firings(2,:);
labels=firings(3,:);

Na=size(pre2a,2);
Nb=size(pre2b,2);

indsA=find(times<Na);
indsB=find(times>=Na);

firingsA=firings(:,indsA);
firingsB=firings(:,indsB); firingsB(2,:)=firingsB(2,:)-Na;

mv.raw=[mfile_path,'/output_s1/pre2.mda'];
mv.firings=firingsA;
mv.sampling_freq=30000;
ms_mountainview(mv);


mv.raw=[mfile_path,'/output_r1/pre2.mda'];
mv.firings=firingsB;
mv.sampling_freq=30000;
ms_mountainview(mv);

end