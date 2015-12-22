function test1

global AA;

X=readmda('example_data/clips_filt2_white.mda');
XI=readmda('example_data/clips_filt2_white_index.mda');
XI=[XI,size(X,3)];
TT=readmda('example_data/templates0_filt2_white.mda');
NN=size(TT,3);
AA=zeros(NN,NN);
for k1=1:NN
for k2=1:NN
    tmp1=TT(:,:,k1); tmp1=tmp1(:);
    tmp2=TT(:,:,k2); tmp2=tmp2(:);
    val=(tmp1'*tmp2)/sqrt(tmp1'*tmp1*tmp2'*tmp2);
    if (val>0.5)
        fprintf('%d,%d\n',k1,k2);
        auc=do_test1(X,XI,k1,k2,1);
        AA(k1,k2)=1-auc;
    end;
end;
end;
figure; imagesc(AA); colormap('gray'); colorbar;

end

function ret=do_test1(X,XI,k1,k2,verbose)

if (k1==k2) ret=0; return; end;

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

mean1=mean(tmp1); stdev1=sqrt(var(tmp1));
mean2=mean(tmp2); stdev2=sqrt(var(tmp2));

stdev0=sqrt((stdev1^2+stdev2^2)/2);
tmp_all=[tmp1+stdev0;tmp2-stdev0];
[~,~,~,auc]=perfcurve([ones(size(tmp1));ones(size(tmp2))*2],tmp_all,2);
%ret=q2-q1;
ret=auc;

if verbose
figure;
plot(tmp1,rand(size(tmp1)),'b.'); hold on;
plot(tmp2,rand(size(tmp2)),'r.'); hold on;
plot([mean1,mean1],ylim,'b'); hold on;
plot([mean1-stdev1,mean1-stdev1],ylim,'r'); hold on;
plot([mean1+stdev1,mean1+stdev1],ylim,'r'); hold on;
plot([mean2,mean2],ylim,'b'); hold on;
plot([mean2-stdev2,mean2-stdev2],ylim,'r'); hold on;
plot([mean2+stdev2,mean2+stdev2],ylim,'r'); hold on;
title(sprintf('k1=%d, k2=%d, auc=%g',k1,k2,auc));
drawnow; pause(0.01);
end;


end
