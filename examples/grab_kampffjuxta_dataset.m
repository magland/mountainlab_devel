function d = grab_kampffjuxta_dataset(n)
% gets Neto-Kampff '16 juxtacellular validated MEA datasets (already in MDA)
%
% d = grab_kampffjuxta_dataset(n)
%  n=1:   2014_11_25_Pair_3_0        (JC is close; neuron easy to get accurate)
%  n=2:   2014_03_26_Pair_2_0        (JC is further away; 5 min)
% Note: these are the two datasets analysed by SpikeDeteck in Neto's paper.
%
% Barnett 4/27/16

if nargin==0, n=1; end

% keep output dir in a fast-IO drive, distinct by expt:
d.outdir = [tempdir,sprintf('kampff%doutput',n)];
if ~exist(d.outdir,'dir'), mkdir(d.outdir); end

mfile_path=fileparts(mfilename('fullpath'));
ext = [mfile_path,'/../ext_datasets'];
dir = [ext,'/Kampff'];

if n==1
  exptname = '2014_11_25_Pair_3_0';
else
  exptname = '2014_03_26_Pair_2_0';
end
dir = [dir '/' exptname];
d.timeseries = [dir,'/raw.mda'];
d.name = ['Kampff-',exptname];
d.samplerate=3e4;
d.truefirings = [dir,'/truefirings.mda'];
d.elecadjmat = [dir,'/elecadjacencymatrix.mda'];    % file, not the array

if 0 % plotting to see if any correlation of Y clips with JC firing:
  Y = readmda(d.timeseries); f = readmda(d.truefirings);
  spikespy({Y,f(2,:),f(3,:),'kampff2'});
end
