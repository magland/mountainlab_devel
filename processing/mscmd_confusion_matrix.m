function mscmd_confusion_matrix(clusters1_path,clusters2_path,output_path,max_matching_offset)

if (nargin==0) test_mscmd_confusion_matrix; return; end;

if (nargin<3) max_matching_offset=3; end;

cmd=sprintf('%s confusion_matrix --clusters1=%s --clusters2=%s --output=%s --max_matching_offset=%d ',mscmd_exe,clusters1_path,clusters2_path,output_path,max_matching_offset);

fprintf('\n*** CONFUSION MATRIX ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end

function test_mscmd_confusion_matrix
T1=[10,20,   30,40,50,60,80   ];
L1=[1, 1,    1, 2, 2, 3, 4    ];
T2=[11,22,21,30,32,52,60      ];
L2=[2, 2, 1, 4, 2  1, 3       ];
C1=zeros(3,length(T1)); C1(2,:)=T1; C1(3,:)=L1;
C2=zeros(3,length(T2)); C2(2,:)=T2; C2(3,:)=L2;
writemda(C1,'tmp_C1.mda')
writemda(C2,'tmp_C2.mda')
opts=struct;
%opts.map12=[1,2,0];
mscmd_confusion_matrix('tmp_C1.mda','tmp_C2.mda','tmp_confusion_matrix.mda',3)
cm=readmda('tmp_confusion_matrix.mda');
disp(cm);
delete('tmp_C1.mda');
delete('tmp_C2.mda');
delete('tmp_confusion_matrix.mda');
end