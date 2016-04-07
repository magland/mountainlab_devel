function d = grab_EJnbhd_dataset()
% returns struct with raw filenames for EJ 7-electrode neighborhood data.
% You need to set up ../ext_datasets as softlink to top of your raw data
% location.
% Barnett 3/21/16

%dir = '~ahb/ss_datasets/EJ/RawMEA2005';
mfile_path=fileparts(mfilename('fullpath'));
dir = [mfile_path,'/../ext_datasets/EJ/RawMEA2005'];
rawfile = '2005-04-26_elec359.mat';
d.outdir = '/tmp/output';
if ~exist(d.outdir,'dir'), mkdir(d.outdir); end
fname = strcat(dir,'/',rawfile);
load(fname);  % is a .mat
d.timeseries = [d.outdir,'/EJnbhd_raw.mda'];
writemda(data, d.timeseries,'float32'); 
d.samplerate = samplingRate;
d.name = 'EJ nbhd M=7 2005-04-26 elec359';
