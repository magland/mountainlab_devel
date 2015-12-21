function test1

global AA;

X=readmda('example_data/clips_filt2_white.mda');
XI=readmda('example_data/clips_filt2_white_index.mda');

do_test1(X,XI,30,31,1); return;

NN=38;

close all;
fA=figure;
AA=zeros(NN,NN);
for k1=1:NN
for k2=1:NN
    fprintf('%d,%d\n',k1,k2);
    AA(k1,k2)=do_test1(X,XI,k1,k2,0);
    figure(fA);
    imagesc(AA); colormap('gray'); colorbar;
    drawnow;
end;
end;




end

function ret=do_test1(X,XI,k1,k2,verbose)
clips1=X(:,:,XI(k1)+1:XI(k1+1));
clips2=X(:,:,XI(k2)+1:XI(k2+1));

[~,~,NC1]=size(clips1);
[~,~,NC2]=size(clips2);

mean1=mean(clips1,3);
mean2=mean(clips2,3);
diff0=mean2-mean1;
diff0=diff0/sum(diff0(:).^2);

tmp1=squeeze(sum(sum(repmat(diff0,1,1,NC1).*clips1,1),2));
tmp2=squeeze(sum(sum(repmat(diff0,1,1,NC2).*clips2,1),2));

median1=median(tmp1);
q1=median1+(median1-median(tmp1(find(tmp1<=median1))));
median2=median(tmp2);
q2=median2-(median(tmp2(find(tmp2>=median2)))-median2);

if verbose
figure;
plot(tmp1,rand(size(tmp1)),'b.'); hold on;
plot(tmp2,rand(size(tmp2)),'r.'); hold on;
plot([q1,q1],ylim,'b'); hold on;
plot([q2,q2],ylim,'r');
title(sprintf('k1=%d, k2=%d',k1,k2));
drawnow; pause(0.01);
end;

ret=q2-q1;

end
