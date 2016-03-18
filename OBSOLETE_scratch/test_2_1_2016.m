function test_2_1_2016

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0=[mfile_path,'/../franklab/test_hippocampal_01_28_2016/tetrode2_output'];

o_detect.threshold=4;
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=30;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=200;

mscmd_detect([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[mfile_path,'/tmp_detect.mda'],[mfile_path,'/tmp_clips.mda'],o_extract_clips);

fprintf('Reading...\n');
detect=readmda([mfile_path,'/tmp_detect.mda']);
clips=readmda([mfile_path,'/tmp_clips.mda']);
pre2=readmda([path0,'/pre2.mda']);
%templates=readmda([path0,'/templates.mda']);
%templates=templates.*abs(templates);

fprintf('Features...\n');
FF=ms_event_features(clips,6);
FF_norms=sqrt(sum(FF.^2,1));
aa=repmat(FF_norms,size(FF,1),1);
FF=FF./aa.*aa.^(1/4);
fprintf('isosplit...\n');
labels=isosplit(FF);
figure; ms_view_clusters(FF,labels);
figure; ms_view_templates_from_clips(clips,labels);
return;

%figure; ms_view_templates(templates);
k=2;
template_k=clips(:,:,randi(NC));
template_k_norm=sqrt(sum(sum(template_k.^2,1),2));

figure; ms_view_templates(template_k);

aaa=clips.*repmat(template_k,1,1,NC);
ips=reshape(sum(sum(aaa,1),2),1,NC);
clip_norms=reshape(sqrt(sum(sum(clips.^2,1),2)),1,NC);

q1=(ips/template_k_norm)/template_k_norm;
q2=sqrt(clip_norms.^2-ips.^2/template_k_norm^2)/template_k_norm;
figure; plot(q1,q2,'b.');


end