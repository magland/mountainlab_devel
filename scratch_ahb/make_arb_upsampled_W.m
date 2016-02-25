% make very upsampled EJ waveforms. Needs validspike tools in path
% Barnett 2/15/16

W = readmda('unit_tests/demo_data/waveforms_EJ_2005-04-26_elec359_60kHz_K7.mda');
[M T K] = size(W);  % check alignment...
tcen = floor((T+1)/2); plot_spike_shapes(W,'W',[],(tcen-1)*ones(1,K));
fac = 10;     % multiplies the 3x upsampling in the above W
[Wu t] = upsample(W,fac);
[M T K] = size(Wu);  % check alignment...
tcen = floor((T+1)/2); plot_spike_shapes(Wu,'Wu',[],(tcen-1)*ones(1,K));
writemda(Wu,'unit_tests/demo_data/waveforms_EJ_2005-04-26_elec359_600kHz_K7.mda');
