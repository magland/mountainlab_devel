function mscmd_confusion_matrix(clusters1_path,clusters2_path,output_path,max_matching_offset)

if (nargin<3) max_matching_offset=3; end;

cmd=sprintf('%s confusion_matrix --clusters1=%s --clusters2=%s --output=%s --max_matching_offset=%d ',mscmd_exe,clusters1_path,clusters2_path,output_path,max_matching_offset);

fprintf('\n*** CONFUSION MATRIX ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end
