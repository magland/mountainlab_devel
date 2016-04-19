% accuracy-measuring function on default dataset, for wrapper to validspike
% Barnett 4/19/16

clear
o = [];
oa.verb=2;  % shows spikespy
[fk Q] = accuracy_anysorter_groundtrutheddata(@validspike_wrapper,[],oa,o);
fprintf('jfm_april_sort sort done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');
fprintf('extended confusion matrix:\n')
for i=1:size(Q,1), fprintf('%d\t',Q(i,:)), fprintf('\n\n'); end
