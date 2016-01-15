function test_detect_001

close all;

mfile_path=fileparts(mfilename('fullpath'));
X=readmda(sprintf('%s/../example_data/filt2_white.mda',mfile_path));

ch=3;
Y=X(ch,:);
if 1
    opts.detect_interval=40;
    opts.detect_threshold=3;
    opts.clip_size=200;
    [Tpos,Tneg]=ms_detect(X(ch,:),opts);
    times=sort([Tpos,Tneg]);
    writemda(times,'scratch/times.mda');
end;
times=readmda('scratch/times.mda');

A=X(:,times);
FF=pca_features(A);
figure; plot(FF(1,:),FF(2,:),'b.'); hold on;

labels=isosplit_mscmd(FF(1:3,:));
max(labels)
ss_view_clusters(FF(1:3,:),labels);

end

function labels=isosplit_mscmd(X)
M=size(X,1);
L=size(X,2);
tmp=zeros(M+2,L);
tmp(1,:)=1;
tmp(2,:)=1;
tmp(3:end,:)=X;
input_path='isosplit_tmp_in.mda';
output_path='isosplit_tmp_out.mda';
writemda(tmp,input_path);
mscmd_cluster(input_path,output_path);
out=readmda(output_path);
delete(input_path);
delete(output_path);
labels=out(3,:);
end

function FF=pca_features(X)
[U,D]=eig(X*X');
[d,I]=sort(diag(D),'descend');
U=U(:,I);
FF=U'*X;
end
