function mscmd_cross_correlograms(clusters_path,output_path,max_dt)

if nargin==0, test_mscmd_cross_correlograms; return; end;

if (nargin<3) max_dt=1500; end;

cmd=sprintf('%s cross_correlograms --clusters=%s --output=%s --max_dt=%d ',mscmd_exe,clusters_path,output_path,max_dt);

fprintf('\n*** CROSS CORRELOGRAMS ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end

function test_mscmd_cross_correlograms

times=[50:50:200,200:20:300];
labels=ones(size(times));
CC=zeros(3,length(times));
CC(2,:)=times;
CC(3,:)=labels;
writemda(CC,'tmp_CC.mda');
mscmd_cross_correlograms('tmp_CC.mda','tmp_cross_correlograms.mda');
cc=readmda('tmp_cross_correlograms.mda');
delete('tmp_CC.mda');
delete('tmp_cross_correlograms.mda');
figure; hist(cc(3,:),1000);

[~,cc2]=ms_cross_correlograms(times,labels,1500);
figure; hist(cc2(3,:),1000);

end