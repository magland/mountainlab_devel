function test_hippocampal_02_18_2016

close all; drawnow;

%%%% Parameters and settings
tetrode_num=1;
shell_opts.min_shell_count=2000;
shell_opts.shell_increment=1;
shell_opts.num_features=12;
shell_opts.merge_threshold=0.8;
o_mild_filter.samplefreq=30000;
o_mild_filter.freq_min=50;
o_mild_filter.freq_max=10000;
o_mild_filter.outlier_threshold=500;
o_filter.samplefreq=30000;
o_filter.freq_min=100;
o_filter.freq_max=10000;
o_filter.outlier_threshold=500;
o_detect.threshold=4;
o_detect.individual_channels=0;
o_detect.normalize=0;
o_detect.inner_window_width=15;
o_detect.outer_window_width=1000;
o_extract_clips.clip_size=120;
o_whiten=struct;

%%%% Set up paths
mfile_path=fileparts(mfilename('fullpath'));
raw_path=[mfile_path,'/../raw/hippocampal/tetrode'];
path0=[mfile_path,sprintf('/output_20151205_tetrode%d',tetrode_num)];
if ~exist(path0,'dir') mkdir(path0); end;

%%%% Extract raw data
extract_raw_data(raw_path,path0);

%%%% Preprocessing
mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre1.mda'],o_filter);
mscmd_bandpass_filter([path0,'/pre0.mda'],[path0,'/pre0_mild.mda'],o_mild_filter);
mscmd_whiten([path0,'/pre1.mda'],[path0,'/pre2.mda'],o_whiten);
mscmd_detect([path0,'/pre2.mda'],[path0,'/detect.mda'],o_detect);
mscmd_extract_clips([path0,'/pre2.mda'],[path0,'/detect.mda'],[path0,'/clips.mda'],o_extract_clips);

%%%% Reading clips
fprintf('Reading...\n');
clips=readmda([path0,'/clips.mda']);
[M,T,NC]=size(clips);
clips=clips-repmat(mean(clips,2),1,T,1); %subtract mean over time

%%%% Shell cluster
fprintf('Shell cluster...\n');
[labels,peaks]=shell_cluster(clips,shell_opts);
K=max(labels);
firings=zeros(4,NC);
detect=readmda([path0,'/detect.mda']);
times=detect(2,:);
firings(1:2,:)=detect;
firings(3,:)=labels;
firings(4,:)=peaks;
writemda(firings,[path0,'/firings.mda']);
%mscmd_templates([path0,'/pre2.mda'],[path0,'/firings.mda'],[path0,'/templates.mda'],o_templates);
%templates=readmda([path0,'/templates.mda']);
%figure; 
%ms_view_templates(templates);

%%%% Writing output and preparing view
fprintf('Writing output and preparing view...\n');
pre2=readmda([path0,'/pre2.mda']);
[clips1,clips1_index]=ms_create_clips_index(ms_extract_clips(pre2,times,o_extract_clips.clip_size),labels);
writemda(clips1,[path0,'/clips0.mda']);
writemda(clips1_index,[path0,'/clips0_index.mda']);
writemda(firings,[path0,'/firings.mda']);
%writemda(corr_matrix,[path0,'/correlation_matrix.mda']);

%%%% Cross correlograms and templates
% mscmd_cross_correlograms([path0,'/firings.mda'],[path0,'/cross_correlograms.mda'],cross_correlograms_max_dt);
% mscmd_templates([path0,'/pre0_mild.mda'],[path0,'/firings.mda'],[path0,'/templates_raw.mda'],struct('clip_size',200));
% mscmd_templates([path0,'/pre2.mda'],[path0,'/firings.mda'],[path0,'/templates.mda'],struct('clip_size',200));
% templates=readmda([path0,'/templates.mda']);
% figure; ms_view_templates(templates);

%%%% MountainView
% view_params.raw=[path0,'/pre2.mda'];
% view_params.firings=[path0,'/firings.mda'];
% view_params.cross_correlograms=[path0,'/cross_correlograms.mda'];
% view_params.templates=[path0,'/templates.mda'];
% view_params.clips=[path0,'/clips0.mda'];
% view_params.clips_index=[path0,'/clips0_index.mda'];
% ms_mountainview(view_params);

view_params.mode='overview2';
view_params.raw=[path0,'/pre2.mda'];
view_params.firings=[path0,'/firings.mda'];
view_params.sampling_freq=o_filter.samplefreq;
ms_mountainview(view_params);

%%%% Split clusters by peak amplitudes
% labels_split=split_clusters_by_peak_amplitudes(clips,labels);
% K_split=max(labels_split);
% templates_split=zeros(M,T,K_split);
% for k=1:K_split
%     inds_k=find(labels_split==k);
%     templates_split(:,:,k)=compute_clips_template(clips(:,:,inds_k));
%     inds0=find_spikes_that_do_not_fit_well(clips(:,:,inds_k),templates_split(:,:,k));
%     fprintf('k=%d: Using %d/%d spikes (%d%%).\n',k,length(inds_k)-length(inds0),length(inds_k),floor((length(inds_k)-length(inds0))/length(inds_k)*100));
%     labels_split(inds_k(inds0))=0;
% end;
% firings_split=readmda([path0,'/firings.mda']);
% firings_split(3,:)=labels_split;
% writemda(firings_split,[path0,'/firings_split.mda']);
% figure; ms_view_templates(templates_split);
% writemda(templates_split,[path0,'/templates_split.mda']);
% mscmd_create_clips_file([path0,'/pre2.mda'],[path0,'/firings_split.mda'],[path0,'/clips0_split.mda'],[path0,'/clips0_split_index.mda'],struct('clip_size',o_extract_clips.clip_size));
% mscmd_cross_correlograms([path0,'/firings_split.mda'],[path0,'/cross_correlograms_split.mda'],cross_correlograms_max_dt);
% 
% %%%% MountainView
% view_params.raw=[path0,'/pre2.mda'];
% view_params.firings=[path0,'/firings_split.mda'];
% view_params.cross_correlograms=[path0,'/cross_correlograms_split.mda'];
% view_params.templates=[path0,'/templates_split.mda'];
% view_params.clips=[path0,'/clips0_split.mda'];
% view_params.clips_index=[path0,'/clips0_split_index.mda'];
% ms_mountainview(view_params);

end

function inds0=find_spikes_that_do_not_fit_well(clips,template)
[M,T,L]=size(clips);
weights=get_template_weights(template,5);
Vclips=reshape(clips,M*T,L);
Vtemplate=reshape(template,M*T,1);
Vweights=reshape(weights,M*T,1);
Vclips=Vclips.*repmat(Vweights,1,L);
Vtemplate=Vtemplate.*Vweights;
Vresid=Vclips-repmat(Vtemplate,1,L);
Vresid_norm=sqrt(sum(Vresid.^2,1));
%Vtemplate_norm=sqrt(Vtemplate'*Vtemplate);
Vweights_norm=sqrt(Vweights'*Vweights);

sigma=1;
expected_vresid_norm=Vweights_norm*sigma;
score=Vresid_norm./expected_vresid_norm;
inds0=find(score>1.5);
end

function Y=get_template_weights(template,num_pix)
[M,T]=size(template);

aa=ifftshift(-floor(T/2):-floor(T/2)+T-1);
sig=num_pix;
kernel=exp(-0.5*aa.^2/sig^2);

fhat=fft(abs(template),[],2);
fhat=fhat.*repmat(kernel,M,1);
Y=real(ifft(fhat,[],2));
end


function template=compute_clips_template(clips)
[M,T,NC]=size(clips);
if (length(clips(:))==0)
    template=zeros(M,T);
    return;
end;
[M,T,NC]=size(clips);
num_features=18;
FF=ms_event_features(clips,num_features);
FFmm=ms_geometric_median(FF);
diffs=FF-repmat(FFmm,1,NC);
dists=sqrt(sum(diffs.^2,1));
sorted_dists=sort(dists);
dist_cutoff=sorted_dists(ceil(length(sorted_dists)*0.3));
inds=find(dists<=dist_cutoff);
template=mean(clips(:,:,inds),3);
end

function [peak_mins,peak_maxs]=define_shells(clip_peaks,opts)

% do the negatives
peak_mins_neg=[];
peak_maxs_neg=[];
if (min(clip_peaks)<0)
    max0=0;
    min0=max0-opts.shell_increment;
    while 1
        if (max0<min(clip_peaks)) break; end;
        if (length(find((min0<=clip_peaks)&(clip_peaks<max0)))>=opts.min_shell_count)
            if (length(find(clip_peaks<min0))>=opts.min_shell_count)
                peak_mins_neg=[peak_mins_neg,min0];
                peak_maxs_neg=[peak_maxs_neg,max0];
                max0=min0;
                min0=max0-opts.shell_increment;
            else
                peak_mins_neg=[peak_mins_neg,-inf];
                peak_maxs_neg=[peak_maxs_neg,max0];
                max0=-inf; min0=-inf;
            end;
        else
            min0=min0-opts.shell_increment;
            if (min0<min(clip_peaks))
                peak_mins_neg=[peak_mins_neg,-inf];
                peak_maxs_neg=[peak_maxs_neg,max0];
                max0=-inf; min0=-inf;
            end;
        end;
    end;
end;
% do the overlap
for j=1:length(peak_maxs_neg)-1
    peak_mins_neg(j)=peak_mins_neg(j+1);
end;

% do the positives
peak_mins_pos=[];
peak_maxs_pos=[];
if (max(clip_peaks)>0)
    min0=0;
    max0=min0+opts.shell_increment;
    while 1
        if (min0>max(clip_peaks)) break; end;
        if (length(find((min0<=clip_peaks)&(clip_peaks<max0)))>=opts.min_shell_count)
            if (length(find(clip_peaks>=max0))>=opts.min_shell_count)
                peak_mins_pos=[peak_mins_pos,min0];
                peak_maxs_pos=[peak_maxs_pos,max0];
                min0=max0;
                max0=min0+opts.shell_increment;
            else
                peak_mins_pos=[peak_mins_pos,min0];
                peak_maxs_pos=[peak_maxs_pos,inf];
                min0=inf; max0=inf;
            end;
        else
            max0=max0+opts.shell_increment;
            if (max0>max(clip_peaks))
                peak_mins_pos=[peak_mins_pos,min0];
                peak_maxs_pos=[peak_maxs_pos,inf];
                min0=inf; max0=inf;
            end;
        end;
    end;
end;
% do the overlap
for j=1:length(peak_mins_pos)-1
    peak_maxs_pos(j)=peak_maxs_pos(j+1);
end;

% combine and sort
peak_mins=[peak_mins_neg,peak_mins_pos];
peak_maxs=[peak_maxs_neg,peak_maxs_pos];
[~,inds]=sort(peak_mins);
peak_mins=peak_mins(inds);
peak_maxs=peak_maxs(inds);

end

function [labels,clip_peaks]=shell_cluster(clips,opts)
[M,T,NC]=size(clips);

fprintf('Computing peaks...\n');
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

[peak_mins,peak_maxs]=define_shells(clip_peaks,opts);

clusterings={};
for ii=1:length(peak_mins)
    rr1=peak_mins(ii);
    rr2=peak_maxs(ii);
    fprintf('rr1=%g rr2=%g... ',rr1,rr2);
    inds_shell=find((clip_peaks>=rr1)&(clip_peaks<rr2));
    CC.inds=inds_shell;
    CC.rr1=rr1;
    CC.rr2=rr2;
    if (length(inds_shell)>1)
        clips_shell=clips(:,:,inds_shell);
        fprintf('features... ');
        FF_shell=ms_event_features(clips_shell,opts.num_features);
        fprintf('isosplit... ');
        labels_shell=isosplit2(FF_shell,struct('whiten_at_each_comparison',1));
        K=max(labels_shell);
        fprintf('K=%d\n',K);
        CC.labels=labels_shell;
        CC.K=K;
    else
        CC.labels=[];
        CC.K=0;
    end;
    clusterings{end+1}=CC;
end

for ii=1:length(clusterings)
    K1=clusterings{ii}.K;
    clusterings{ii}.clusters={};
    for k1=1:K1
        clusterings{ii}.clusters{k1}.inds=clusterings{ii}.inds(find(clusterings{ii}.labels==k1));
    end;
    if ii>1
        K2=clusterings{ii-1}.K;
        match_counts=zeros(K1,K2);
        for k1=1:K1
            for k2=1:K2
                inds_k1=clusterings{ii}.inds(find(clusterings{ii}.labels==k1));
                inds_k2=clusterings{ii-1}.inds(find(clusterings{ii-1}.labels==k2));
                match_counts(k1,k2)=length(intersect(inds_k1,inds_k2));
            end;
        end;
        for k1=1:K1
            for k2=1:K2
                numer=match_counts(k1,k2);
                denom=sum(match_counts(k1,:))+sum(match_counts(:,k2))-match_counts(k1,k2);
                if ((denom)&&(numer/denom>=opts.merge_threshold))
                    pct=numer/denom;
                    fprintf('Merging [%d,%d] to [%d,%d] (%d%%)\n',ii,k1,ii-1,k2,floor(pct*100));
                    inds_k1=clusterings{ii}.clusters{k1}.inds;
                    inds_k2=clusterings{ii-1}.clusters{k2}.inds;
                    clusterings{ii}.clusters{k1}.inds=union(inds_k1,inds_k2);
                    clusterings{ii-1}.clusters{k2}.inds=[];
                end;
            end;
        end;
    end;
end;

clusters={};
for ii=1:length(clusterings)
    for jj=1:length(clusterings{ii}.clusters)
        CC=clusterings{ii}.clusters{jj};
        if (length(CC.inds)>0)
            clusters{end+1}=CC;
        end;
    end;
end;

labels=zeros(1,NC);
for jj=1:length(clusters)
    labels(clusters{jj}.inds)=jj;
end;

K=max(labels);
used=zeros(1,K); used(labels)=1; iii=find(used);
mapping=zeros(1,K);
for jj=1:length(iii) mapping(iii(jj))=jj; end;
labels=mapping(labels);
K=max(labels);

end

function extract_raw_data(raw_path,output_path)

raw_mda_fname=sprintf('%s/20151205_s1_tet1_refchan.mda',raw_path);
tetrode_fname=sprintf('%s/pre0.mda',output_path);

if (~exist(tetrode_fname,'file'))
    fprintf('Reading raw data...\n');
    raw=readmda(raw_mda_fname);
    raw=raw';

    tetrode=raw(1:4,:)-repmat(raw(5,:),4,1);
    fprintf('Writing tetrode data...\n');
    writemda(tetrode,tetrode_fname);
end;

L=[...
    0,4;...
    0,3;...
    0,2;...
    0,1;...
];
writemda(L,[output_path,'/locations.mda']);

end

function [labels_out]=split_clusters_by_peak_amplitudes(clips,labels)
[M,T,NC]=size(clips);
templates=zeros(M,T,0);
K=max(labels);
labels_out=zeros(1,NC);
for k=1:K
    inds=find(labels==k);
    labels0=split_cluster_by_peak_amplitudes(clips(:,:,inds));
    labels_out(inds)=labels0+max(labels_out);
end;
end

function [labels]=split_cluster_by_peak_amplitudes(clips)
[M,T,NC]=size(clips);
clip_peaks_pos=squeeze(max(clips(:,T/2+1,:),[],1))';
clip_peaks_neg=-squeeze(max(-clips(:,T/2+1,:),[],1))';
clip_peaks=clip_peaks_pos.*(abs(clip_peaks_pos)>abs(clip_peaks_neg))+clip_peaks_neg.*(abs(clip_peaks_pos)<abs(clip_peaks_neg));

[peak_mins,peak_maxs]=define_shells(clip_peaks,struct('shell_increment',2,'min_shell_count',200));

templates=zeros(M,T,length(peak_mins));
labels=zeros(1,NC);
for ii=1:length(peak_mins)
    inds=find((clip_peaks>=peak_mins(ii))&(clip_peaks<peak_maxs(ii)));
    labels(inds)=max(labels)+1;
end;

end
