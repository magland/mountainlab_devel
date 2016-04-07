function d = grab_harris2000_dataset()
% Get classic d5331 dataset w/ IC channel output as true firings. Barnett 4/6/16

mfile_path=fileparts(mfilename('fullpath'));
ext = [mfile_path,'/../ext_datasets'];
dir = [ext,'/Harris2000/d5331'];
rawfile = 'd533101.dat';
fname = strcat(dir,'/',rawfile);
fname
fid = fopen(fname,'r');
if fid==-1, error('Harris data not found!'); end
Y = fread(fid,inf,'int16');
fclose(fid);
n = numel(Y)
Nch = 8;
Y = reshape(Y,[Nch n/Nch]);
imagesc(Y);

% *** pull out the IC



Y = Y(2:5,:);  % the EC channels
d.outdir = '/tmp/output'; if ~exist(d.outdir,'dir'), mkdir(d.outdir); end  % fix
d.timeseries = [d.outdir,'/harris2000_raw.mda'];
writemda(Y, d.timeseries,'float32'); 
d.name = 'Harris2000 d5331';
d.samplerate = 1e4;
