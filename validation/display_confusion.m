function display_confusion(Q,aname,bname)
% DISPLAY_CONFUSION   pretty plots and analysis of confusion matrix
%
%  display_confusion(Q,aname,bname) makes big plot of confusion matrix Q, ro- and
%   col-normalized, and recall/precision descending plots. aname and bname
%   are the (optional) names for the horizontal and vertical axes.
%
% Barnett 5/13/16

if nargin<2, aname = 'A'; end
if nargin<3, bname = 'B'; end
na = sum(Q,2); nb = sum(Q,1); Ka = numel(na); Kb = numel(nb); K = min(Ka,Kb);
mKb = size(Q,2)-1;
Qa = Q./repmat(na,[1 numel(nb)]); % row sum normalize
Qb = Q./repmat(nb,[numel(na) 1]); % col "
figure; set(gcf,'position',[1100 400 3000 1000]);  % huge, 4k monitor!
subplot(1,6,1:2);
imagesc(Qa); colormap(1-gray(256)); title('best Q, row-normalized');
ylabel([aname ' label']); xlabel([bname ' label']);
hold on; plot([.5,mKb+.5;mKb+1.5,mKb+.5], [Ka+.5,.5;Ka+.5,Ka+1.5],'r-');
subplot(1,6,3:4);
imagesc(Qb); colormap(1-gray(256)); title('best Q, col-normalized');
ylabel([aname ' label']); xlabel([bname ' label']);
hold on; plot([.5,mKb+.5;mKb+1.5,mKb+.5], [Ka+.5,.5;Ka+.5,Ka+1.5],'r-');
subplot(1,6,5);
recall = diag(Q)./na(1:K); [srecall,sri] = sort(recall,'descend'); % sri=sorted indices
plot(srecall,'+-'); text(1:K,srecall,num2cellstr(sri)); title('sorted recalls');
axis([1 K 0 1]); xlabel('# of labels');
subplot(1,6,6);
prec = diag(Q)./nb(1:K)'; [sprec,spi] = sort(prec,'descend');
plot(sprec,'+-');text(1:K,sprec,num2cellstr(spi)); title('sorted precisions');
axis([1 K 0 1]); xlabel('# of labels');
drawnow
