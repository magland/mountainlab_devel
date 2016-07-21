% analysis of sorting Kampff datasets. Assumes sorting already done.
% truth from juxtacellular elec.
% Barnett 7/21/16

% To sort:  /mnt/xfs1/home/ahb/ss_datasets/Kampff/alg_scda_005_kampff.sh

clear
% the easy one...
dir = '/mnt/xfs1/home/ahb/ss_datasets/Kampff/2014_11_25_Pair_3_0';
d.timeseries = [dir,'/pre2.mda'];
d.firings = [dir,'/firings.mda'];
d.truefirings = [dir,'/truefirings.mda'];
d.name = 'kampff 2014_11_25_Pair_3_0';
d.samplerate = 3e4;    % from READ ME.docx
o_acc.verb = 2;
[fk Q perm] = accuracy_groundtruthedsorting(d,o_acc)

% thresh    |   acc fk                     | events 
%   3.0          0.994  (3 mislabeled, 1 falspos, out of 348)     750000
% This compares vs Neto '16 p.10:  38 false pos, no missed.


if 1   % the hard one...
dir = '/mnt/xfs1/home/ahb/ss_datasets/Kampff/2014_03_26_Pair_2_0';
d.timeseries = [dir,'/pre2.mda'];
d.firings = [dir,'/firings.mda'];
d.truefirings = [dir,'/truefirings.mda'];
d.name = 'kampff 2014_03_26_Pair_2_0';
d.samplerate = 3e4;    % from READ ME.docx
o_acc.verb = 2;
[fk Q perm] = accuracy_groundtruthedsorting(d,o_acc)
end

% finds 58 out of 150.  (beats Neto which finds 35 / 150)
% but 5549 false pos !


mv.mode='overview2'; mv.raw=d.timeseries; mv.samplerate=3e4;
mv.firings = d.firings;  mountainview(mv);
