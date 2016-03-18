% Show how to call the simple sorter on the demo data.
% Run from the mountainsort top directory.
% Barnett 3/18/16

clear
forceregen = 0;
[Yfile,~,~,o.samplerate] = get_default_dataset([],forceregen); % Yfile: EC signal
output_dir = 'unit_tests/demo_data/output';
if ~exist(output_dir,'dir'), mkdir(output_dir); end

[firingsfile,prefile]=simplesorter(Yfile,output_dir,o);

% load and view the EC input with firings output...
Y = readmda(Yfile);
firings = readmda(firingsfile); times=firings(2,:); labels=firings(3,:);
spikespy({Y,times,labels,'simple sorter'});
K = max(labels); pops = histc(labels,1:K); disp('populations n_k vs k:');
fprintf('\t%d',1:K); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');
