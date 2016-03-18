function test_filter

close all;

mfile_path=fileparts(mfilename('fullpath'));
path0=sprintf('%s/../compare_ms_kk/output_ms',mfile_path);

templates=readmda([path0,'/templates0_filt2_white.mda']);

ss_view_waveforms(templates);

[M,T,K]=size(templates);


ch=7;

NN=1e6;
YY=zeros(1,NN);
for j=1:1500
    kk=randi(K);
    tt=randi(NN-T);
    YY(1,tt:tt+T-1)=YY(1,tt:tt+T-1)+templates(ch,:,kk);
end;

YYhat=fft(YY);
figure; plot(abs(YYhat(1:end)));

end