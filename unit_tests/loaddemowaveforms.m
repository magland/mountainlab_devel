function [W samplerate] = loaddemowaveforms
% [W samplerate] = loaddemowaveforms returns MxTxK templates array in W,
%    and sample rate in Hz in samplerate.

% Barnett 2/19/16. Upsampled 2/24/16

% not upsampled case:
%demowaveforms = 'unit_tests/demo_data/waveforms_EJ_2005-04-26_elec359_20kHz_K7.mda';
%samplerate = 2e4;

% 3x upsampled:
%demowaveforms = 'unit_tests/demo_data/waveforms_EJ_2005-04-26_elec359_60kHz_K7.mda';
%samplerate = 6e4;

% 30x upsampled:
demowaveforms = 'unit_tests/demo_data/waveforms_EJ_2005-04-26_elec359_600kHz_K7.mda';
samplerate = 6e5;

W = readmda(demowaveforms);
