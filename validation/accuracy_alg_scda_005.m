% accuracy-measuring function on default dataset, for jfm alg_scda_005
% Barnett 7/28/16

clear
o.sign = -1;      % sorting opts
o.detectmeth=3;  % 1 or 3

oa.verb=1;  % 2 shows spikespy
[fk Q perm info] = accuracy_anysorter_groundtrutheddata(@alg_scda_005,[],oa,o);
fprintf('alg_scda_005 sort done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');
fprintf('extended confusion matrix:\n')
for i=1:size(Q,1), fprintf('%d\t',Q(i,:)), fprintf('\n'); end
