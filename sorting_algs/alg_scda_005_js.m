function [firingsfile,info] = alg_scda_005_js(tsfile,path,o)
% ALG_SCDA_005_JS  Matlab wrapper to JFM's javascript/C++ sorter of July 2017
%
% [firingsfile,info] = alg_scda_005_js(rawfile,output_dir,o)
%
% Inputs:
%    rawfile - path to .mda of MxN (ie # channels by # timepoints) raw signal data
%    output_dir - path to existing directory where all output will be written
%    o - sorting options:
%            o.parfile: path of JSON parameter file (overrides all other options,
%                       required)
%            o.eleccoordsfile: path of CSV file giving x,y coords of each electrode
%                       (optional; if absent assumes full dense connectivity)
% Outputs:
%    firingsfile - path to the firings.mda output file
%    info - struct with fields:
%           filtfile - filtered timeseries
%           prefile - path to the preprocessed timeseries (filt and whitened)
%
% Also see: ml/matlab/sorting_algorithms/alg_scda_005.m  which is MATLAB
%           reimplementation of this.

mfile_path=fileparts(mfilename('fullpath'));
ml = [mfile_path, '/../../mountainlab/'];
if isfield(o,'eleccoordsfile')  % specify electrode geom...
  cmd = sprintf('%s/mountainprocess/bin/mountainprocess run-script --_nodaemon %s/sorting_algorithms/alg_scda_005.js %s --raw=%s --geom=%s --outpath=%s',...
                ml,ml,o.parfile,tsfile,o.eleccoordsfile,path)
else       % assumed full connectivity graph for electrodes...
  cmd = sprintf('%s/mountainprocess/bin/mountainprocess run-script --_nodaemon %s/sorting_algorithms/alg_scda_005.js %s --raw=%s --outpath=%s',...
                ml,ml,o.parfile,tsfile,path)
end
system(cmd);
firingsfile = [path,'/firings.mda'];
info.filtfile = [path,'/pre1b.mda'];
info.prefile = [path,'/pre2.mda'];
