function mscmd_fit(X_path,cluster_in_path,templates_path,cluster_out_path)

if (nargin<1) example_mscmd; return; end;

X=readmda(X_path);
cluster=readmda(cluster_in_path);
templates=readmda(templates_path);

%sort by time
[~,sort_inds]=sort(cluster(2,:));
cluster=cluster(:,sort_inds);

num_events=size(cluster,2);
clip_size=size(templates,2);
tt1=-floor(clip_size/2);
tt2=tt1+clip_size-1;

to_include=zeros(1,num_events);
times=cluster(2,:);
labels=cluster(3,:);

while 1
    fprintf('pass... ');
    scores=compute_scores(X,cluster,templates);
    num_changed=0;
    for j=1:num_events
        if ((~to_include(j))&&(scores(j)>0))
            best=1;
            k=j-1;
            while (k>=1)&&(abs(times(k)-times(j))<=clip_size)
                if (scores(k)>scores(j)) best=0; break; end;
                k=k-1;
            end;
            if (best)
                k=j+1;
                while (k<=num_events)&&(abs(times(k)-times(j))<=clip_size)
                    if (scores(k)>scores(j)) best=0; break; end;
                    k=k+1;
                end;
            end;
            if (best)
                if ((times(j)+tt1>=1)&&(times(j)+tt2<=size(X,2)))
                    to_include(j)=1;
                    num_changed=num_changed+1;
                    X(:,times(j)+tt1:times(j)+tt2)=X(:,times(j)+tt1:times(j)+tt2)-templates(:,:,labels(j));
                end;
            end;
        end;
    end;
    fprintf('%d events added\n',num_changed);
    if (num_changed==0) break; end;
end;

inds=find(to_include);
cluster_out=cluster(:,inds);

writemda(cluster_out,cluster_out_path);

end

function scores=compute_scores(X,cluster,templates)

num_events=size(cluster,2);
clip_size=size(templates,2);
tt1=-floor(clip_size/2);
tt2=tt1+clip_size-1;

times=cluster(2,:);
labels=cluster(3,:);

scores=zeros(1,num_events);

for j=1:num_events
if ((times(j)+tt1>=1)&&(times(j)+tt2<=size(X,2)))
    X0=X(:,times(j)+tt1:times(j)+tt2);
    template0=templates(:,:,labels(j));
    tmp1=X0;
    tmp2=X0-template0;
    sumsqr1=sum(tmp1(:).^2);
    sumsqr2=sum(tmp2(:).^2);
    scores(j)=sumsqr1-sumsqr2;
end;
end;

end
