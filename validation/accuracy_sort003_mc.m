% Showing how to use the accuracy-measuring function on demo sorter 001
% Barnett 3/17/16

clear; addpath('franklab/sort_003_multichannel');
o = [];
%o.detect_threshold = 3.0;   % ds001 uses thresh in sigma units
fk = accuracy_anysorter_groundtrutheddata(@sort_003_multichannel,[],'/tmp/output',[],o);
fprintf('Franklab 003 MC sort done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');

