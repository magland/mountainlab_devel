function view_waveforms_spikespy(waveforms)
[M,T,K]=size(waveforms);
padding=ceil(T/2);
times0=(0:K-1)*(T+padding);
labels0=1:K;
waveforms0=reshape(cat(2,waveforms,zeros(M,padding,K)),[M,(T+padding)*K]);
spikespy({waveforms0,times0,labels0});
end
