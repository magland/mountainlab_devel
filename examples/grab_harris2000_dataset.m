function d = grab_harris2000_dataset(raw_fname)
% Get classic d5331 dataset w/ IC channel output as true firings.
% Barnett 4/6/16
% Changed excluded interval 8/9/16

mfile_path=fileparts(mfilename('fullpath'));

if nargin<1
    ext = [mfile_path,'/../ext_datasets'];
    dir = [ext,'/Harris2000/d5331'];
    rawfile = 'd533101.dat';
    fname = strcat(dir,'/',rawfile);
else
    fname=raw_fname;
end;

d.outdir = '/tmp/harrisoutput'; if ~exist(d.outdir,'dir'), mkdir(d.outdir); end  % fix

fprintf('reading %s ...\n',fname)
fid = fopen(fname,'r');
if fid==-1, error('Harris data not found!'); end
Y = fread(fid,inf,'int16');
fclose(fid);
n = numel(Y);
Nch = 8;
N = n/Nch; if N~=round(N), error('non-integer number of time pts for given M'); end
Y = reshape(Y,[Nch N]);
d.samplerate = 1e4;
t = (1:N)/d.samplerate;   % time indices
%j = find((t>26 & t<109.49) | t>110.29);   % ahb original cutout excluded opening
j = find(t<109.49 | t>110.29);   % cut out 0.8 sec of bad bursting/noise part, forgets absolute t. As in Ekanadham et al 2013.
N = numel(j);
Y = Y(:,j);

YIC = Y(6,:); % pull out the IC (intra-cellular), for ground truth
%writemda(YIC-mean(YIC),[d.outdir,'/harris2000_IC.mda'],'float32'); % for debug
trig = (max(YIC)+min(YIC))/2;   % trigger level (checked by eye)
truetimes = 1+find(diff(YIC>trig)==1); % upwards-going times, sample units
%figure; plot(1:N, YEC); vline(truetimes); xlabel('t (in samples)');
truefirings = [0*truetimes; truetimes; 1+0*truetimes];  % chan; time; label
d.truefirings = [d.outdir,'/harris2000_truefirings.mda'];
writemda(truefirings, d.truefirings,'float64');

Y = Y(2:5,:);  % the EC channels
d.timeseries = [d.outdir,'/harris2000_raw.mda'];
writemda(Y, d.timeseries,'float32'); 
d.dims = size(Y);
d.name = 'Harris2000 d5331';
