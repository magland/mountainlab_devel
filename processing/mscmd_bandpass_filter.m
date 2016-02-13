function mscmd_bandpass_filter(input_path,output_path,opts)
%MSCMD_BANDPASS_FILTER - Bandpass filter using smooth roll-offs in Fourier space
%
%The is the mountainsort version of ms_filter
%
% Syntax:  mscmd_bandpass_filter(input_path,output_path,opts)
%
% Inputs:
%    input_path - mda file of MxN input raw data
%    output_path - mda file of MxN output filtered data
%    opts.samplefreq - the sampling frequency corresponding to X, e.g.
%                      30000
%    opts.freq_min - the lower end of the bandpass filter
%    opts.freq_max - the upper end of the bandpass filter
%
%
% Other m-files required: mscmd_exe
%
% See also: ms_filter, spikespy

% Author: Jeremy Magland
% Jan 2016; Last revision: 13-Feb-2016

if (~isfield(opts,'outlier_threshold')) opts.outlier_threshold=0; end;

cmd=sprintf('%s bandpass_filter --input=%s --output=%s ',mscmd_exe,input_path,output_path);
cmd=[cmd,sprintf('--samplefreq=%g --freq_min=%g --freq_max=%g --outlier_threshold=%g',opts.samplefreq,opts.freq_min,opts.freq_max,opts.outlier_threshold)];

fprintf('\n*** BANDPASS FILTER ***\n');
fprintf('%s\n',cmd);
status=system(cmd);

if (status~=0)
    error('mountainsort returned with error status %d',status);
end;

end