function mscmd_cross_correlograms(firings_path,output_path,max_dt)
%MSCMD_CROSS_CORRELOGRAMS - Generate cross-correlogram data, which in
%particular may be used for MountainView visualization.
%
%MountainView can be used to view cross-correlograms. It requires the
%cross-correlograms to be prepared with a particular format. This function
%outputs a file that is prepared to be read by MountainView
%
%This is a fast implementation and is much preferred over
%ms_cross_correlograms.
%
% Syntax:  mscmd_cross_correlograms(firings_path,output_path,max_dt)
%
% Inputs:
%    firings_path - path of input array of times/labels. See docs for the
%                    format of the firings.mda file.
%    output_path - path of the output .mda file that will contain the array
%                  of cross-correlogram data, to be read in by MountainView
%    max_dt - the maximum integer timepoint interval to consider
%
%
% Other m-files required: mscmd_exe
%
% See also: ms_cross_correlograms, ms_mountainview

% Author: Jeremy Magland
% Jan 2016; Last revision: 13-Feb-2016

if nargin==0, test_mscmd_cross_correlograms; return; end;

if (nargin<3) max_dt=1500; end;

cmd=sprintf('%s cross_correlograms --firings=%s --output=%s --max_dt=%d ',mscmd_exe,firings_path,output_path,max_dt);

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