function mscmd_branch_cluster_v1(raw_path,detect_path,adjacency_matrix_path,out_firings_path,opts)

if (nargin==0) test_mscmd_branch_cluster_v1; return; end;

if (nargin<5) opts=struct; end;

if (~isfield(opts,'clip_size')) opts.clip_size=100; end;
if (~isfield(opts,'min_section_count')) opts.min_section_count=50; end;
if (~isfield(opts,'section_increment')) opts.section_increment=0.5; end;
if (~isfield(opts,'num_features')) opts.num_features=3; end;

cmd=sprintf('%s branch_cluster_v1 --raw=%s --detect=%s --adjacency_matrix=%s --firings=%s --clip_size=%d --min_section_count=%d --section_increment=%g --num_features=%d',mscmd_exe,raw_path,detect_path,adjacency_matrix_path,out_firings_path,...
    opts.clip_size,opts.min_section_count,opts.section_increment,opts.num_features);

fprintf('\n*** BRANCH_CLUSTER_V1 ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end

function test_mscmd_branch_cluster_v1

close all;

rng(1);
raw_path=default_test_dataset('synthetic1');

if ~exist('tmp_output','dir') mkdir('tmp_output'); end;
[firings_path,pre_path]=sort_002_multichannel(raw_path,'tmp_output');
mv.raw=pre_path;
mv.firings=firings_path;
mv.sampling_freq=20000;
ms_mountainview(mv);

firings=readmda(firings_path);
channels=firings(1,:);
times=firings(2,:);
labels=firings(3,:);
peaks=firings(4,:);

clips=ms_extract_clips(readmda(pre_path),times,100);
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

mscmd_bandpass_filter(raw_path,'tmp_pre1.mda',struct('samplefreq',20000,'freq_min',100,'freq_max',10000));
mscmd_whiten('tmp_pre1.mda','tmp_pre2.mda');
mscmd_detect('tmp_pre2.mda','tmp_detect.mda',struct('detect_threshold',4,'detect_interval',50,'individual_channels',1,'clip_size',100));
detect=readmda('tmp_detect.mda');
mscmd_branch_cluster_v1('tmp_pre2.mda','tmp_detect.mda','','tmp_firings.mda');
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

mv.raw='tmp_pre2.mda';
mv.firings='tmp_firings.mda';
mv.sampling_freq=20000;
ms_mountainview(mv);

return;

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