
% AHB: this seems to be all funcs for JFM's self-contained MATLAB clustering
% as of 2/10/16

function [labels,clip_peaks]=jfm_shell_cluster(clips,opts)
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

%%%%%%%%%%%%%%%%%%%%%%%%


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
