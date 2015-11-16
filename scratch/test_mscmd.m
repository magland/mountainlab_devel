addpath('mountainview/src/spikespy/matlab');

path0='example1_data';
opts.samplefreq=30000;
opts.freq_min=300;
opts.freq_max=10000;
opts.ncomp=4;

oo.num_channels=72; oo.channels=[37:52,68,69];
oo.t1=0; oo.t2=19e6; oo.threshold=500;
mscmd_extract([path0,'/ms11d45.dat'],[path0,'/test0.mda'],oo);

mscmd_bandpass_filter([path0,'/test0.mda'],[path0,'/test1.mda'],opts);
%mscmd_normalize_channels([path0,'/test1.mda'],[path0,'/test2.mda']);
mscmd_whiten([path0,'/test1.mda'],[path0,'/test3.mda'],opts);

X=readmda_data_beginning([path0,'/test1.mda'],3e6);
spikespy(X);

X=readmda_data_beginning([path0,'/test3.mda'],3e6);
spikespy(X);