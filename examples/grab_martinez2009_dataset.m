function d = grab_martinez2009_dataset(n)
% Get Martinez-Quiroga 2009 JNM 1-chan sim datasets w/ known firings.
%
% d = grab_martinez2009_dataset(n) where n = 1...5 returns struct pointing
% to MDA conversion from one of 5 datasets.

% Barnett 4/11/16

d.outdir = '/tmp/output'; if ~exist(d.outdir,'dir'), mkdir(d.outdir); end  % fix

mfile_path=fileparts(mfilename('fullpath'));
ext = [mfile_path,'/../ext_datasets'];
dir = [ext,'/MartinezQuiroga2009_sims'];
rawfile = sprintf('simulation-%d');
fname = strcat(dir,'/',rawfile);
fprintf('reading %s ...\n',fname)
load(fname);



d.samplerate = 1e4;
t = (1:N)/d.samplerate;   % time indices
j = find((t>26 & t<109.49) | t>110.29);   % cut out bad parts & forget absolute t
N = numel(j);
Y = Y(:,j);

YEC = Y(6,:); % pull out the IC (intra-cellular), for ground truth
trig = (max(YEC)+min(YEC))/2;   % trigger level (checked by eye)
truetimes = 1+find(diff(YEC>trig)==1); % upwards-going times, sample units
%figure; plot(1:N, YEC); vline(truetimes); xlabel('t (in samples)');
truefirings = [0*truetimes; truetimes; 1+0*truetimes];  % chan; time; label
d.truefirings = [d.outdir,'/harris2000_truefirings.mda'];
writemda(truefirings, d.truefirings,'float64');

Y = Y(2:5,:);  % the EC channels
d.timeseries = [d.outdir,'/harris2000_raw.mda'];
writemda(Y, d.timeseries,'float32'); 
d.name = 'Harris2000 d5331';

