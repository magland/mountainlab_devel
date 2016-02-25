function [W samplerate] = loaddemowaveforms
% [W samplerate] = loaddemowaveforms returns MxTxK templates array in W,
%    and sample rate in Hz in samplerate.

% Barnett 2/19/16. Upsampled 2/24/16

% not upsampled case:
demowaveforms = 'unit_tests/demo_data/waveforms_EJ_2005-04-26_elec359_20kHz_K7.mda';
W = readmda(demowaveforms);
samplerate = 2e4;

% upsampled:
demowaveforms = 'unit_tests/demo_data/waveforms_EJ_2005-04-26_elec359_60kHz_K7.mda';
W = readmda(demowaveforms);
samplerate = 6e4;
