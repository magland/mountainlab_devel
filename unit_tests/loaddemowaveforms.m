function [W samplerate] = loaddemowaveforms
% [W samplerate] = loaddemowaveforms
%
% Barnett 2/19/16
demowaveforms = 'unit_tests/demo_data/waveforms_EJ_2005-04-26_elec359_20kHz_K7.mda';
W = readmda(demowaveforms);
samplerate = 2e4;
