% Show how to call the simple sorter on the demo data.
% Barnett 3/18/16

clear
d = demo_dataset;   % get struct pointing to demo data files

opts.samplerate = d.samplerate;
[firingsfile,~] = simplesorter(d.signal,d.outdir,opts);

% load and view the EC input signal with firings output...
Y = readmda(d.signal);
firings = readmda(firingsfile); times=firings(2,:); labels=firings(3,:);
spikespy({Y,times,labels,'simple sorter'});
K = max(labels); pops = histc(labels,1:K); disp('populations n_k vs k:');
fprintf('\t%d',1:K); fprintf('\n'); fprintf('\t%d',pops); fprintf('\n');
