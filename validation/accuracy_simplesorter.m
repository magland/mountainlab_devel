% Showing how to use the accuracy-measuring function on the simple sorter
% Barnett 3/16/16

forceregen = 1;
d = demo_dataset(forceregen);
o.detect_threshold = 4.0;   % sorter opts: threshold in stddev EJ demo dataset
o.detect_polarity='m';
fk = accuracy_anysorter_groundtrutheddata(@simplesorter,d,[],o);
fprintf('accuracy of simplesorter done; fk accuracies vs label k are... \n')
fprintf('\t%d',1:numel(fk)), fprintf('\n'); fprintf('\t%.3f',fk), fprintf('\n');

