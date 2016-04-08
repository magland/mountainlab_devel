% accuracy-measuring function on default dataset, for jfm april sort.
% Barnett 4/8/16

clear
o.sign = -1;
[fk Q] = accuracy_anysorter_groundtrutheddata(@jfm_april_sort,[],[],o);
fprintf('jfm_april_sort sort done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');
fprintf('extended confusion matrix:\n')
for i=1:size(Q,1), fprintf('%d\t',Q(i,:)), fprintf('\n\n'); end
