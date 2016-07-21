function [h Qa Qb] = display_confusion(Q,aname,bname)
% DISPLAY_CONFUSION   pretty plots and analysis of confusion matrix
%
%  display_confusion(Q,aname,bname) makes big plot of the extended confusion
%   matrix Q, row- and col-normalized, and recall & precisions in descending
%   orders (with labels added as text numbers). The latter two plots only
%   are meaningful when Q is input in best-permuted ordering.
%   aname and bname are the (optional) names for the vertical and horizontal
%   axes respectively
%
%  h = display_confusion(...) outputs array of handles to axes

% Barnett 5/13/16; tweaks 6/8/16

if nargin<2, aname = 'A'; end
if nargin<3, bname = 'B'; end
na = sum(Q,2); nb = sum(Q,1); Ka = size(Q,1)-1; Kb = size(Q,2)-1; K = min(Ka,Kb);
Qa = Q./repmat(na,[1 numel(nb)]); % row sum normalize
Qb = Q./repmat(nb,[numel(na) 1]); % col "

figure; set(gcf,'position',[1100 400 3000 1000]);  % huge, 4k monitor!

h(1) = subplot(1,6,1:2);
imagesc(Qa); colormap(1-gray(256)); hold on;
popsc = 2e3/(max(Ka,Kb)*sqrt(max(na)));   % blob size scale (propto area)
for k=1:Ka+1
  if na(k)>0, plot(0,k,'b.','markersize',max(1,popsc*sqrt(na(k)))); end
end
axis([-.5 Kb+1.5 -.5 Ka+1.5]);
set(gca,'ytick',1:Ka+1)
set(gca,'yticklabel',[num2cellstr(1:Ka); '\phi'])
title('extended confusion matrix, row-normalized');
ylabel([aname ' label']); xlabel([bname ' label']);
hline(Ka+.5,'r-'); vline(Kb+.5,'r-'); hline(.5,'k-'); vline(.5,'k-');

h(2) = subplot(1,6,3:4);
imagesc(Qb); colormap(1-gray(256)); hold on
popsc = 2e3/(max(Ka,Kb)*sqrt(max(nb)));   % blob size scale (propto area)
for k=1:Kb+1
  if nb(k)>0, plot(k,0,'b.','markersize',max(1,popsc*sqrt(nb(k)))); end
end
axis([-.5 Kb+1.5 -.5 Ka+1.5]);
set(gca,'xtick',1:Kb+1)
set(gca,'xticklabel',[num2cellstr(1:Kb); '\phi'])
title('extended confusion matrix, col-normalized');
ylabel([aname ' label']); xlabel([bname ' label']);
hline(Ka+.5,'r-'); vline(Kb+.5,'r-'); hline(.5,'k-'); vline(.5,'k-');

h(3) = subplot(1,6,5);
recall = diag(Q(1:K,1:K))./na(1:K);
recall(isnan(recall)) = 0;            % kill nans so they move to end
[srecall,sri] = sort(recall,'descend'); % sri=sorted indices
plot(srecall,'.-'); text(0.2+(1:K),0.01+srecall,num2cellstr(sri));
title('recalls (descending order)');
axis([1 max(2,K) 0 1]); xlabel('# of matching labels');

h(4) = subplot(1,6,6);
prec = diag(Q(1:K,1:K))./nb(1:K)';
prec(isnan(prec)) = 0;            % kill nans so they move to end
[sprec,spi] = sort(prec,'descend');
plot(sprec,'.-');text(0.2+(1:K),0.01+sprec,num2cellstr(spi));
title('precisions (descending order)');
axis([1 max(2,K) 0 1]); xlabel('# of matching labels');

drawnow

