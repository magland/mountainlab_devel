function ms_view_templates_from_clips(clips,clips_index_or_labels,opts)

if nargin<1, test_ms_view_templates_from_clips; return; end;

if (nargin<3) opts=struct; end;
if ~isfield(opts,'clip_size') opts.clip_size=0; end;
if ~isfield(opts,'show_stdev') opts.show_stdev=0; end;


[M,T,NC]=size(clips);
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

templates=zeros(M,T,K);
if opts.show_stdev
    opts.stdev=zeros(M,T,K);
end;

for k=1:K
    fprintf('%d ',k); if (mod(k,10)==0) fprintf('\n'); end;
    ii0=find(labels==k);
    clips0=clips(:,:,ii0);
    templates(:,:,k)=mean(clips0,3);
    if opts.show_stdev
        opts.stdev(:,:,k)=sqrt(var(clips0,[],3));
    end;
end;
fprintf('\n');

ms_view_templates(templates,opts);

end

function test_ms_view_templates_from_clips

mfile_path=fileparts(mfilename('fullpath'));
path0=sprintf('%s/../example_data',mfile_path);

clips=readmda([path0,'/clips_filt.mda']);
cii=readmda([path0,'/clips_filt_index.mda']);

opts.show_stdev=1;

figure;
ms_view_templates_from_clips(clips,cii,opts);

end
