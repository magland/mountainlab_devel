function ms_view_clip_clouds(clips,clips_index_or_labels,opts)

if nargin<1, test_ms_view_clip_clouds; return; end;

[M,T,NC]=size(clips);

if (nargin<3) opts=struct; end;
if ~isfield(opts,'clip_size') opts.clip_size=0; end;
if ~isfield(opts,'vzoom') opts.vzoom=1; end;
if (~isfield(opts,'Tpad')) opts.Tpad=floor(T*0.2); end;

Tpad=opts.Tpad;

if (opts.clip_size>0)
    clips=clips(:,T/2-opts.clip_size/2:T/2+opts.clip_size/2-1,:);
    [M,T,NC]=size(clips);
end;

if (length(clips_index_or_labels)==NC)
    labels=clips_index_or_labels;
else
    clips_index=clips_index_or_labels;
    K=length(clips_index);
    clips_index=[clips_index,NC];
    labels=[];
    for k=1:K
        ii0=(clips_index(k)+1):clips_index(k+1);
        labels=[labels,ones(1,length(ii0))*k];
    end;
end;

K=max(labels);

nbins=600;
tmp=max(clips,[],2);
stdev0=sqrt(var(tmp(:)));
ybins=linspace(-stdev0*4,stdev0*4,nbins);
ybins=ybins/opts.vzoom; %zoom vertically
XX=zeros(T+Tpad,K,length(ybins),M);

for k=1:K
    fprintf('%d ',k); if (mod(k,10)==0) fprintf('\n'); end;
    ii0=find(labels==k);
    clips0=clips(:,:,ii0);
    for tt=1:T
        for mm=1:M
            tmp=hist(squeeze(clips0(mm,tt,:)),ybins)/length(ii0);
            XX(tt,k,:,mm)=tmp(end:-1:1);
        end;
    end;
end;
fprintf('\n');
XX=reshape(XX,(T+Tpad)*K,length(ybins)*M);

imagesc(XX');

set(gca,'xtick',[]); set(gca,'ytick',[]);
set(gca,'position',[0,0,1,1]);
colormap('parula');
set(gca,'clim',get(gca,'clim')*0.1); %make it brighter!

end

function test_ms_view_clip_clouds

mfile_path=fileparts(mfilename('fullpath'));
path0=sprintf('%s/../example_data',mfile_path);

clips=readmda([path0,'/clips_filt.mda']);
cii=readmda([path0,'/clips_filt_index.mda']);

opts=struct;

figure;
ms_view_clip_clouds(clips,cii,opts);

end
