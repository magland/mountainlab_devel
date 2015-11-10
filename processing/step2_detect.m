function step2_detect(opts,data)

timerA=tic;

detect_times_prefix=[opts.working_path,'/detect/times_'];
detect_clips_prefix=[opts.working_path,'/detect/clips_'];
adjacency_path=[opts.output_path,'/adjacency.mda'];

X=data.X;

AM=readmda(adjacency_path);

fprintf('Detecting...\n');
oo.thresh=opts.detection_threshold;
oo.interval=opts.detection_interval;
oo.clip_size=opts.clip_size;
for j=1:size(X,1)
    fname_pos=[detect_times_prefix,sprintf('pos_%d.mda',j)];
    fname_neg=[detect_times_prefix,sprintf('neg_%d.mda',j)];
    fname2_pos=[detect_clips_prefix,sprintf('pos_%d.mda',j)];
    fname2_neg=[detect_clips_prefix,sprintf('neg_%d.mda',j)];
    if ((~exist(fname_pos,'file'))||(~exist(fname_neg,'file'))||(~exist(fname2_pos,'file'))||(~exist(fname2_neg,'file')))
        
        Xfilt=ss_freqfilter(X(j,:),30000,opts.detect_freq_min,opts.detect_freq_max);
        for j=1:size(Xfilt,1)
            stdev0=sqrt(var(Xfilt(j,:)));
            Xfilt(j,:)=Xfilt(j,:)/stdev0;
        end;
        
        [Tpos,Tneg]=detect_events(Xfilt,oo);
        
        clips_pos=extract_clips(X,Tpos,opts.clip_size);
        clips_neg=extract_clips(X,Tneg,opts.clip_size);
        
        AM(j,j)=0;
        neighbor_channels=[find(AM(:,j))];
        
%         clips_pos_filt=pca_filter_clips(clips_pos(neighbor_channels,:,:),20);
%         clips_neg_filt=pca_filter_clips(clips_neg(neighbor_channels,:,:),20);
%         
%         inds0=find_valid_clips(clips_pos_filt,opts.detection_threshold,1); Tpos=Tpos(inds0); clips_pos=clips_pos(:,:,inds0);
%         inds0=find_valid_clips(clips_neg_filt,opts.detection_threshold,1); Tneg=Tneg(inds0); clips_neg=clips_neg(:,:,inds0);
        
        fprintf('Writing %s and %s (%d pos, %d neg events detected)...\n',fname_pos,fname_neg,length(Tpos),length(Tneg));
        writemda(Tpos,fname_pos);
        writemda(Tneg,fname_neg);
        
        fprintf('Writing %s and %s...\n',fname2_pos,fname2_neg);
        writemda(clips_pos,fname2_pos);
        writemda(clips_neg,fname2_neg);
    else
        fprintf('Files exists %s, %s, %s and %s...\n',fname_pos,fname_neg,fname2_pos,fname2_neg);
    end;
end;

fprintf('\nElapsed: %g seconds',toc(timerA));
fprintf('\n');

end

function clips_filt=pca_filter_clips(clips,npca)

[M,T,NC]=size(clips);

[FF,info]=ss_eventfeatures(clips);
FF=FF(1:npca,:);
SS=info.subspace;
SS=SS(:,:,1:npca);

SS0=reshape(SS,[M*T,size(SS,3)]);
clips_filt=reshape(SS0*FF,[M,T,NC]);

end

function inds=find_valid_clips(clips,thresh,num)

% CC=(abs(clips)>=thresh)*1.0;
% counts=squeeze(sum(sum(CC,1),2));
% inds=find(counts>=num);

clips=reshape(clips,[size(clips,1)*size(clips,2),size(clips,3)]);
clips=sort(abs(clips),1);
vals=mean(clips(end-num+1:end,:),1);
inds=find(vals>=thresh);

end

function clips=extract_clips(Y,times,clip_size)
% Extract clips centered at times

[M,N]=size(Y);
T=clip_size;
C=length(times);

clips=zeros(M,T,C);
tt1=-ceil((clip_size)/2);
tt2=tt1+clip_size-1;
if (min(times+tt1)<1) error('Invalid time in extract_clips'); end;
if (max(times+tt2)>N) error('Invalid time in extract_clips'); end;
for j=1:C
	clips(:,:,j)=Y(:,times(j)+tt1:times(j)+tt2);
end;

end

function [Tpos,Tneg]=detect_events(X,opts)

interval=opts.interval;
thresh=opts.thresh;

absX=abs(X);
absX=max(absX,1); %max over channels (but best to put in one channel at a time)
times=find(absX>thresh);
signed_vals=X(times);
pos_inds=find(signed_vals>=0);
neg_inds=find(signed_vals<0);

%Now make sure we only use those that are global maxima over the radius
%Nt/2
use_it=ones(1,length(times));
for j=1:length(times)
    if ((times(j)-interval)<1) use_it(j)=0; end;
    if ((times(j)+interval)>size(X,2)) use_it(j)=0; end;
    k=j-1;
    while (k>=1)&&(times(k)>=times(j)-interval)
        if (absX(times(k))>absX(times(j)))
            use_it(j)=0;
        else
            use_it(k)=0;
        end;
        k=k-1;
    end
end;

times_pos=times(pos_inds);
times_neg=times(neg_inds);
Tpos=times_pos(find(use_it(pos_inds)));
Tneg=times_neg(find(use_it(neg_inds)));

% clips_pos=extract_clips(X,Tpos,opts.clip_size);
% stdevs00=squeeze(sqrt(sum(sum(clips_pos.^2,1),2)/(size(clips_pos,1)*size(clips_pos,2))));
% inds00=find(absX(Tpos)>thresh*stdevs00');
% Tpos=Tpos(inds00);
% clips_neg=extract_clips(X,Tneg,opts.clip_size);
% stdevs00=squeeze(sqrt(sum(sum(clips_neg.^2,1),2)/(size(clips_neg,1)*size(clips_neg,2))));
% inds00=find(absX(Tneg)>thresh*stdevs00');
% Tneg=Tneg(inds00);

end

