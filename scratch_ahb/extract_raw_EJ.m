function info = extract_raw_EJ(head,indir,outdir)
% EXTRACT_RAW_EJ.  Converts EJ raw data to MDA raw file.
%
% info = extract_raw_EJ(head,indir,outdir)
% Writes pre0, samplerate, locations
%
% Barnett 2/16/16
load([indir '/' head '.mat']);
fprintf('raw data size: %d by %d\n',size(data,1),size(data,2))
writemda(data,[outdir '/pre0.mda']);  % JFM filenames
writemda(samplingRate,[outdir '/samplerate.mda']);

z = [0, exp(2i*pi*(1:6)/6)]; electrodelocs = [real(z);imag(z)];
writemda(electrodelocs',[outdir '/locations.mda']);    % MDA is M-by-2
info = [];
