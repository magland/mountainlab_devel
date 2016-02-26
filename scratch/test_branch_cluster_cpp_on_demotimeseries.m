function test_branch_cluster_cpp_on_demotimeseries

close all;

%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'));
addpath([mfile_path,'/../franklab/sort_002_multichannel']);
raw_path=[mfile_path,'/../unit_tests/demo_data/demotimeseries.mda'];
if ~exist(raw_path), writedemotimeseries; end

mfile_path=fileparts(mfilename('fullpath'));
output_path=[mfile_path,'/output_sort001_on_demotimeseries'];
if ~exist(output_path,'dir') mkdir(output_path); end;

%%%% Sort
sort_opts=struct;
sort_opts.detectability_threshold=2.5; %this controls exclusion of noise clusters
%[firings_path,pre_path]=sort_002_multichannel(raw_path,output_path,sort_opts);
[firings_path,pre_path]=do_the_sorting(raw_path,output_path,sort_opts);

%%%% View output
mv.mode='overview2';
mv.raw=pre_path;
mv.firings=firings_path;
ms_mountainview(mv);

end

function [firings_path,pre_path]=do_the_sorting(raw_path,output_path,sort_opts)

mscmd_bandpass_filter(raw_path,'tmp_pre1.mda',struct('samplefreq',20000,'freq_min',100,'freq_max',10000));
mscmd_whiten('tmp_pre1.mda','tmp_pre2.mda');
mscmd_detect('tmp_pre2.mda','tmp_detect.mda',struct('detect_threshold',4,'detect_interval',50,'individual_channels',1,'clip_size',100));
mscmd_branch_cluster_v1('tmp_pre2.mda','tmp_detect.mda','','tmp_firings.mda');

detect=readmda('tmp_detect.mda');
channels=detect(1,:);
times=detect(2,:);
pre2=readmda('tmp_pre2.mda');
clips=ms_extract_clips(pre2,times,100);
inds_2=find(channels==2);
clips_2=clips(:,:,inds_2);
FF_2=ms_event_features(clips_2,3);
labels_2=isosplit2(FF_2);
figure; ms_view_clusters(FF_2,labels_2);


firings_path='tmp_firings.mda';
pre_path='tmp_pre2.mda';

firings=readmda('tmp_firings.mda');
channels=firings(1,:);
times=firings(2,:);
labels=firings(3,:);
peaks=firings(4,:);

clips=ms_extract_clips(readmda('tmp_pre2.mda'),times,100);
[M,T,L]=size(clips);
for m=1:M
    inds=find(channels==m);
    if (length(inds)>0)
        FF=ms_event_features(clips(:,:,inds),3);
        LL=labels(inds); LL(LL~=0)=LL(LL~=0)-min(LL(LL~=0))+1;
        figure; ms_view_clusters(FF,LL);
        title(sprintf('m=%d',m));
        figure; ms_view_templates_from_clips(clips(:,:,inds),LL);
        title(sprintf('m=%d',m));
    end;
end;

return;

mv.raw='tmp_pre2.mda';
mv.firings='tmp_firings.mda';
mv.sampling_freq=20000;
ms_mountainview(mv);


end