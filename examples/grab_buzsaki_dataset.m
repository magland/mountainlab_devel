function d = grab_buzsaki_dataset(len)
% d = grab_buzsaki_dataset(len)
% returns struct with raw filenames for Buzsaki 10-channel probe.
% You need to set up ../ext_datasets as softlink to top of your raw data
% location.
% If len absent, uses default
% len = 1: 50 sec,   len = 2: 150 sec,   len = 3: 1500 sec (25 min)
% Barnett 3/24/16

if nargin<1, len = 1; end
mfile_path=fileparts(mfilename('fullpath'));
dir = [mfile_path,'/../ext_datasets/Buzsaki'];
if len==1
  d.name = 'Buzsaki CingulateCortex BXRat19 gp5 unfilt N=1e6';
  head = 'CingulateCortex_BXRat19_gp5_unfilt_1e6';
elseif len==2
  d.name = 'Buzsaki CingulateCortex BXRat19 gp5 unfilt N=3e6';
  head = 'CingulateCortex_BXRat19_gp5_unfilt_3e6';
elseif len==3
  d.name = 'Buzsaki CingulateCortex BXRat19 gp5 unfilt N=3e7';
  head = 'CingulateCortex_BXRat19_gp5_unfilt_3e7';
else
  error('unknown len (try 1,2,3)');
end
d.outdir = '/tmp/output';
%d.outdir = [dir,'/output'];   % is on ceph so is super slow to seek
if ~exist(d.outdir,'dir'), mkdir(d.outdir); end
fname = strcat(dir,'/',head,'.dat');
d.timeseries = [d.outdir,'/',head,'.mda'];
if ~exist(d.timeseries,'file')
  fid=fopen(fname);
  raw = fread(fid,[10 inf],'float'); fclose(fid); % read as many rows as exist
  writemda(raw, d.timeseries,'float32');
end
d.samplerate = 2e4;

