function test_sort002_on_demotimeseries

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
sort_opts.test_mode=0;
[firings_path,pre_path]=sort_002_multichannel(raw_path,output_path,sort_opts);

%%%% View output
mv.mode='overview2';
mv.raw=pre_path;
mv.firings=firings_path;
mv.sampling_freq=20000;
ms_mountainview(mv);

% firings=readmda(firings_path);
% times=firings(2,:);
% labels=firings(3,:);
% raw=readmda(pre_path);
% clips=ms_extract_clips(raw,times,100);
% templates=ms_templates(clips,labels);
% figure; ms_view_templates(templates);
% for m=1:size(templates,1)
%     figure;
%     for k=1:size(templates,3)
%         plot(templates(m,:,k)); hold on;
%     end;
%     title(sprintf('channel %d',m));
% end;

end
